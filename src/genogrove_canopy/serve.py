# SPDX-License-Identifier: GPL-3.0-or-later
"""``canopy serve`` — a local web front-end over the same generate → execute → render pipeline.

One self-contained HTML page (query box, cohort picker, history, generated-code panel, a
records table, and a locus view for interval answers) talks to two endpoints on a stdlib
``http.server``:

* ``GET  /``        — the page
* ``GET  /cohorts`` — the rE2G cohort list, for the picker
* ``POST /ask``     — ``{question, cohort}`` → ``{summary, records, code}`` (or ``{error, code}``)

No web framework: the backend is glue around :mod:`genogrove_canopy.llm` (code-gen) and a warm
:class:`genogrove_canopy.sandbox.Worker` per cohort selection — the same worker ``--interactive`` uses,
so the grove ``open`` is paid once per cohort and reused across questions.
"""

from __future__ import annotations

import argparse
import json
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from genogrove_canopy import llm, resources, sandbox
from genogrove_canopy.cli import (
    DEFAULT_COHORT,
    DEFAULT_MODEL,
    _BASE,
    _cohort_ids,
    _grove_context,
    _parse_output,
    _pygenogrove_site_dir,
    _resolve_cohorts,
)


class _Grove:
    """The plain-GENCODE grove made runnable: system prompt, code preamble, and warm worker.

    Cohort-independent now — the grove is GENCODE structure only, so there's exactly one, built
    once and reused across questions. Enhancers are resolved per question and injected as
    ``ENHANCERS`` (see :func:`_pipeline`), not baked into the grove.
    """

    def __init__(self, model: str) -> None:
        block, self.preamble, data_paths = _grove_context()
        self.system_prompt = llm.build_system_prompt(block)
        self.model = model
        self.worker = sandbox.Worker(
            data_paths=data_paths, extra_syspath=[_pygenogrove_site_dir()]
        )


# The single plain-GENCODE grove, built lazily and reused. A lock serializes /ask (one worker,
# one query at a time). ponytail: global lock — a local single-user tool.
_GROVE: _Grove | None = None
_LOCK = threading.Lock()

#: Cap on an /ask body — a question, not a payload. Guards against a client-supplied
#: Content-Length that would otherwise be read straight into memory.
_MAX_BODY = 64 * 1024


def _grove(model: str) -> _Grove:
    global _GROVE
    if _GROVE is None:
        _GROVE = _Grove(model)
    return _GROVE


def _serve_cohorts(cohort_override: str, cohort_hint: str):
    """Grounded per-query cohort: UI picker override → the model's declared cohort → default.
    Returns ``({name: accessions}, note_or_None)``; no catalog match → ``({}, note)``."""
    if cohort_override:
        return _resolve_cohorts([cohort_override]), None
    if cohort_hint:
        try:
            return _resolve_cohorts([cohort_hint]), None
        except SystemExit:
            return {}, f"no cohort matched {cohort_hint!r} — no enhancers loaded"
    return _resolve_cohorts([DEFAULT_COHORT]), "default"


def _pipeline(question: str, cohort: str, model: str, emit) -> dict:
    """Run one question end to end, calling ``emit(kind, msg)`` before each visible step. The UI
    ``cohort`` is an **override**; otherwise the model's declared cohort (or the default) is used.
    Returns the result dict (``records``/``summary``/``code``/``note``) or ``{"error", "code"?}``."""
    if not resources._all_grove_gg(_BASE).exists():
        emit("step", "Fetching the grove — first run only: a pinned ~109 MB .gg (genes + cCREs)")
        resources.ensure_all_grove(_BASE)

    emit("step", "Opening the grove")
    grove = _grove(model)

    emit("step", "Generating the pygenogrove query")
    cohort_hint, targets, code = llm.generate_query(question, grove.system_prompt, model=grove.model)

    enh_pre, note = "", None
    if targets:
        from genogrove_canopy.layers import enhancers
        # `why` is the internal reason flag; `note` is display text. Keeping them apart stops the
        # "default" sentinel leaking to the UI when a cohort resolves but finds no links.
        cohorts, why = _serve_cohorts(cohort, cohort_hint)
        note = None if why == "default" else why
        cohort_ids = _cohort_ids(cohorts)
        if cohort_ids:
            emit("step", f"Fetching enhancers for cohort(s) {'; '.join(cohorts)}")
            records = enhancers.fetch_for_targets(targets, cohort_ids)
            if records:
                enh_pre = enhancers.preamble(records)
                note = f"{len(records)} enhancer links from {'; '.join(cohorts)}" + (
                    " (default)" if why == "default" else "")

    emit("step", "Running the query over the grove")
    result = grove.worker.submit("import json\n" + grove.preamble + enh_pre + code)
    if result.returncode != 0 or result.timed_out:
        err = result.stderr.strip() or "(the generated code failed with no output)"
        return {"error": err, "code": code}
    records, passthrough = _parse_output(result.stdout)
    return {"summary": passthrough, "records": records, "code": code, "note": note}


class _Handler(BaseHTTPRequestHandler):
    model = DEFAULT_MODEL

    def log_message(self, *args) -> None:  # quiet the default per-request stderr spam
        pass

    def _send(self, code: int, body: str, ctype: str) -> None:
        payload = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:
        if self.path == "/" or self.path.startswith("/?"):
            self._send(200, PAGE, "text/html; charset=utf-8")
        elif self.path == "/cohorts":
            cohorts = [
                {"name": c["name"], "ontology_id": c["ontology_id"],
                 "type": c["type"], "n_replicates": c["n_replicates"]}
                for c in resources.re2g_cohorts()
            ]
            self._send(200, json.dumps({"cohorts": cohorts, "default": DEFAULT_COHORT}),
                       "application/json")
        else:
            self._send(404, "not found", "text/plain")

    def do_POST(self) -> None:
        if self.path != "/ask":
            self._send(404, "not found", "text/plain")
            return
        # Parse before responding, so a malformed request gets a 400 rather than a dead connection:
        # a bad Content-Length, non-JSON body, or non-string field all raise here, and this runs
        # outside the try below (which can only report errors once the 200 stream has started).
        try:
            length = int(self.headers.get("Content-Length", 0))
            if length > _MAX_BODY:
                self._send(413, json.dumps({"error": "request too large"}), "application/json")
                return
            req = json.loads(self.rfile.read(length) or "{}")
            question, cohort = req.get("question") or "", req.get("cohort") or ""
            if not isinstance(question, str) or not isinstance(cohort, str):
                raise TypeError("question and cohort must be strings")
            question, cohort = question.strip(), cohort.strip()
        except (ValueError, TypeError, AttributeError) as exc:
            self._send(400, json.dumps({"error": f"bad request: {exc}"}), "application/json")
            return
        if not question:
            self._send(400, json.dumps({"error": "empty question"}), "application/json")
            return

        # Stream newline-delimited JSON: {t:"step"|"tick", msg} events as work progresses,
        # then a final {t:"done", result} (or {t:"error", msg}). HTTP/1.0 close-delimited —
        # the client reads the body stream to EOF.
        self.send_response(200)
        self.send_header("Content-Type", "application/x-ndjson")
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()

        def event(obj) -> None:
            try:
                self.wfile.write((json.dumps(obj) + "\n").encode("utf-8"))
                self.wfile.flush()
            except (BrokenPipeError, ConnectionError):
                pass  # client navigated away mid-run; nothing to do

        try:
            with _LOCK:  # one build/query at a time — see the note on _LOCK
                result = _pipeline(question, cohort, self.model, lambda kind, msg: event({"t": kind, "msg": msg}))
            event({"t": "done", "result": result})
        except Exception as exc:  # surface the failure in the stream, never crash the server
            event({"t": "error", "msg": f"{type(exc).__name__}: {exc}"})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="canopy serve", description="Local web front-end for canopy."
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8787)
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"Anthropic model for code generation (default: {DEFAULT_MODEL}).")
    parser.add_argument("--no-browser", action="store_true",
                        help="Do not open a browser window on start.")
    args = parser.parse_args(argv)

    _Handler.model = args.model
    server = ThreadingHTTPServer((args.host, args.port), _Handler)
    url = f"http://{args.host}:{args.port}"
    print(f"canopy serve — {url}  (Ctrl-C to stop)")
    print("The first question of each cohort pays the grove open (and, for a new cohort, its "
          "build); after that queries are warm.")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping…")
    finally:
        server.server_close()
        if _GROVE is not None:
            _GROVE.worker.close()
    return 0


# ---------------------------------------------------------------------------
# The page. One file: markup + styles + the record-driven renderer. Kept inline
# so `canopy serve` ships with no static-asset directory to locate at runtime.
# ---------------------------------------------------------------------------
PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>canopy</title>
<style>
  :root {
    color-scheme: light;
    --surface:#fcfcfb; --panel:#ffffff; --ink:#14130f; --ink2:#52514e; --muted:#898781;
    --grid:#e6e5df; --axis:#c9c8c1; --border:rgba(20,19,15,0.10); --border2:rgba(20,19,15,0.16);
    --gene:#6c6a63; --gene-fill:#d8d6cd; --variant:#d8382f;
    --promoter:#2a78d6; --genic:#eb6834; --intergenic:#1baf7a; --feat:#7a67c9;
    --exon:#3f3d37; --exon-utr:#b7b5ac; --zoom-wash:rgba(216,56,47,0.07);
    --accent:#1baf7a; --accent-ink:#0b5c40;
    --tip-bg:#1f1e1b; --tip-ink:#f6f5f0; --code-bg:#f6f5f1;
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) {
      color-scheme: dark;
      --surface:#161513; --panel:#1d1c19; --ink:#f4f3ee; --ink2:#c3c2b7; --muted:#8f8d85;
      --grid:#2c2b28; --axis:#3b3a36; --border:rgba(255,255,255,0.10); --border2:rgba(255,255,255,0.17);
      --gene:#9a988f; --gene-fill:#34332e; --variant:#f06a5f;
      --promoter:#3987e5; --genic:#e06a3a; --intergenic:#22b183; --feat:#9085e9;
      --exon:#d8d6cd; --exon-utr:#55534c; --zoom-wash:rgba(240,106,95,0.10);
      --accent:#22b183; --accent-ink:#bff0dc;
      --tip-bg:#f4f3ee; --tip-ink:#1d1c19; --code-bg:#211f1b;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --surface:#161513; --panel:#1d1c19; --ink:#f4f3ee; --ink2:#c3c2b7; --muted:#8f8d85;
    --grid:#2c2b28; --axis:#3b3a36; --border:rgba(255,255,255,0.10); --border2:rgba(255,255,255,0.17);
    --gene:#9a988f; --gene-fill:#34332e; --variant:#f06a5f;
    --promoter:#3987e5; --genic:#e06a3a; --intergenic:#22b183; --feat:#9085e9;
    --exon:#d8d6cd; --exon-utr:#55534c; --zoom-wash:rgba(240,106,95,0.10);
    --accent:#22b183; --accent-ink:#bff0dc;
    --tip-bg:#f4f3ee; --tip-ink:#1d1c19; --code-bg:#211f1b;
  }
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; }
  body {
    background:var(--surface); color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI",sans-serif; line-height:1.5;
    -webkit-font-smoothing:antialiased;
    display:grid; grid-template-columns:236px 1fr; grid-template-rows:auto 1fr;
    grid-template-areas:"brand top" "side main"; min-height:100vh;
  }
  .mono { font-family:ui-monospace,"SF Mono",Menlo,monospace; font-variant-numeric:tabular-nums; }

  .brand { grid-area:brand; padding:18px 18px 14px; border-bottom:1px solid var(--border);
           border-right:1px solid var(--border); display:flex; align-items:center; gap:9px; }
  .brand .logo { width:22px; height:22px; }
  .brand b { font-size:18px; font-weight:650; letter-spacing:-0.01em; }
  .brand small { color:var(--muted); font-size:11px; letter-spacing:.12em; text-transform:uppercase; }

  .top { grid-area:top; padding:14px 22px; border-bottom:1px solid var(--border);
         display:flex; gap:10px; align-items:center; flex-wrap:wrap; }
  .top form { display:flex; gap:10px; align-items:center; flex:1; min-width:280px; }
  #q {
    flex:1; min-width:220px; padding:11px 14px; font-size:15px; color:var(--ink);
    background:var(--panel); border:1px solid var(--border2); border-radius:10px;
  }
  #q:focus { outline:2px solid var(--accent); outline-offset:1px; border-color:transparent; }
  select {
    padding:10px 12px; font-size:13px; color:var(--ink); background:var(--panel);
    border:1px solid var(--border2); border-radius:10px; max-width:230px;
  }
  button.ask {
    padding:11px 20px; font-size:15px; font-weight:600; color:#fff; background:var(--accent-ink);
    border:none; border-radius:10px; cursor:pointer;
  }
  button.ask:hover { filter:brightness(1.08); }
  button.ask:disabled { opacity:.55; cursor:progress; }

  .side { grid-area:side; border-right:1px solid var(--border); padding:16px 12px; overflow-y:auto; }
  .side h3 { font-size:11px; letter-spacing:.13em; text-transform:uppercase; color:var(--muted);
             margin:2px 8px 10px; font-weight:600; }
  .hist { list-style:none; margin:0; padding:0; display:flex; flex-direction:column; gap:3px; }
  .hist li { padding:8px 10px; border-radius:8px; font-size:13px; color:var(--ink2); cursor:pointer;
             border:1px solid transparent; }
  .hist li:hover { background:var(--panel); border-color:var(--border); color:var(--ink); }
  .hist .empty { color:var(--muted); font-size:12.5px; padding:8px 10px; cursor:default; }

  main { grid-area:main; padding:26px 30px 80px; overflow-y:auto; max-width:1080px; }
  .placeholder { color:var(--muted); font-size:15px; max-width:60ch; margin-top:8px; }
  .placeholder .ex { display:block; margin-top:12px; }
  .placeholder .ex code { display:inline-block; margin:4px 8px 0 0; padding:5px 10px; font-size:13px;
    background:var(--panel); border:1px solid var(--border); border-radius:8px; cursor:pointer; color:var(--ink2); }
  .placeholder .ex code:hover { color:var(--ink); border-color:var(--border2); }

  .status { font-size:13px; color:var(--muted); margin:0 0 14px; min-height:18px; }
  .status.err { color:var(--variant); }
  .note { color:var(--accent-ink); }

  .qecho { font-size:20px; font-weight:600; letter-spacing:-0.01em; margin:0 0 4px; text-wrap:balance; }
  .summary { background:var(--panel); border:1px solid var(--border); border-radius:10px;
             padding:12px 14px; margin:14px 0; font-size:13.5px; white-space:pre-wrap; }

  .progress { margin-top:18px; }
  .steps { list-style:none; margin:0; padding:0; display:flex; flex-direction:column; gap:9px; }
  .steps li { display:flex; align-items:center; gap:11px; font-size:13.5px; color:var(--ink); }
  .steps li.done { color:var(--muted); }
  .steps li .lbl { flex:0 1 auto; }
  .steps li .tick { color:var(--accent-ink); font-size:12px; font-variant-numeric:tabular-nums; }
  .steps li .spin { width:14px; height:14px; flex:none; border-radius:50%;
    border:2px solid var(--grid); border-top-color:var(--accent); animation:spin .8s linear infinite; }
  .steps li.done .spin { border:none; animation:none; position:relative; }
  .steps li.done .spin::after { content:"✓"; position:absolute; inset:0; text-align:center;
    color:var(--accent); font-size:13px; line-height:14px; }
  @keyframes spin { to { transform:rotate(360deg); } }
  @media (prefers-reduced-motion:reduce){ .steps li .spin { animation:none; } }

  section.block { margin-top:26px; }
  .cap { display:flex; align-items:baseline; justify-content:space-between; gap:12px; margin-bottom:9px; }
  .cap h2 { font-size:14px; font-weight:650; margin:0; }
  .cap p { margin:0; font-size:12px; color:var(--muted); }
  .plot { background:var(--panel); border:1px solid var(--border); border-radius:12px; padding:8px 6px; overflow-x:auto; }
  svg { display:block; width:100%; height:auto; }
  svg text { font-family:ui-monospace,Menlo,monospace; font-variant-numeric:tabular-nums; }

  .legend { display:flex; flex-wrap:wrap; gap:15px; margin-top:12px; font-size:12.5px; color:var(--ink2); }
  .legend span { display:inline-flex; align-items:center; gap:7px; }
  .legend i { width:22px; height:3px; border-radius:2px; display:inline-block; }
  .legend i.round { width:11px; height:11px; border-radius:50%; }
  .legend i.vline { width:3px; height:13px; border-radius:1px; }
  .legend i.blk { width:14px; height:9px; border-radius:2px; }

  details.code { margin-top:22px; border:1px solid var(--border); border-radius:10px; background:var(--panel); }
  details.code summary { cursor:pointer; padding:11px 14px; font-size:13px; font-weight:600; color:var(--ink2); }
  details.code[open] summary { border-bottom:1px solid var(--border); }
  details.code pre { margin:0; padding:14px; overflow-x:auto; font-size:12.5px; line-height:1.55;
    font-family:ui-monospace,Menlo,monospace; background:var(--code-bg); border-radius:0 0 10px 10px; }

  .tablewrap { overflow-x:auto; border:1px solid var(--border); border-radius:10px; margin-top:10px; }
  table { border-collapse:collapse; width:100%; font-size:12.5px; }
  th,td { text-align:left; padding:7px 11px; white-space:nowrap; border-bottom:1px solid var(--border); }
  th { color:var(--muted); font-weight:600; position:sticky; top:0; background:var(--panel);
       font-size:11px; letter-spacing:.04em; text-transform:uppercase; }
  td.mono { font-variant-numeric:tabular-nums; }
  tbody tr:last-child td { border-bottom:none; }
  tbody tr:hover td { background:var(--code-bg); }
  .swatch { display:inline-block; width:8px; height:8px; border-radius:50%; margin-right:6px; vertical-align:baseline; }

  #tip { position:fixed; z-index:30; pointer-events:none; opacity:0; transition:opacity .08s;
    background:var(--tip-bg); color:var(--tip-ink); padding:8px 10px; border-radius:8px;
    font-size:12px; max-width:280px; box-shadow:0 6px 24px rgba(0,0,0,.28); line-height:1.45; }
  #tip .t-h { font-weight:650; margin-bottom:2px; }
  #tip .t-r { font-family:ui-monospace,Menlo,monospace; opacity:.85; font-size:11px; }
  .hoverable { cursor:default; }
  .hoverable:hover { filter:brightness(1.09); }
  @media (prefers-reduced-motion:reduce){ #tip{transition:none;} }
  @media (max-width:720px){
    body { grid-template-columns:1fr; grid-template-areas:"brand" "top" "main"; }
    .side { display:none; }
  }
</style>
</head>
<body>
  <div class="brand">
    <svg class="logo" viewBox="0 0 24 24" fill="none" aria-hidden="true">
      <path d="M12 21V11" stroke="var(--gene)" stroke-width="2" stroke-linecap="round"/>
      <path d="M12 12c-3.5 0-6-2.2-6-5.5C9.5 6.5 12 8.7 12 12z" fill="var(--intergenic)"/>
      <path d="M12 10c0-3.3 2.5-5.8 6-5.8C18 7.5 15.5 10 12 10z" fill="var(--genic)"/>
      <circle cx="12" cy="21" r="1.8" fill="var(--variant)"/>
    </svg>
    <span><b>canopy</b><br><small>locus view</small></span>
  </div>

  <div class="top">
    <form id="askform">
      <input id="q" name="q" autocomplete="off" placeholder="Ask about a variant, gene, or its enhancers…">
      <select id="cohort" title="Enhancer cohort (rE2G)"></select>
      <button class="ask" type="submit">Ask</button>
    </form>
  </div>

  <aside class="side">
    <h3>History</h3>
    <ul class="hist" id="hist"><li class="empty">Your questions appear here.</li></ul>
  </aside>

  <main id="main">
    <p class="qecho" style="font-weight:600">Ask a question over connected genomic intervals.</p>
    <div class="placeholder">
      canopy translates plain English to a <span class="mono">pygenogrove</span> query, runs it over the
      GENCODE grove (plus a chosen enhancer→gene layer), and draws where the answer lands.
      <span class="ex">
        <code>Which gene contains the variant at chr7:55,191,822?</code>
        <code>Which transcripts have an exon over chr7:55,191,822?</code>
        <code>What are the enhancers of EGFR?</code>
      </span>
    </div>
  </main>

  <div id="tip"></div>

<script>
const SVGNS = "http://www.w3.org/2000/svg";
const $ = s => document.querySelector(s);
const CLS = { promoter:"var(--promoter)", genic:"var(--genic)", intergenic:"var(--intergenic)" };

function el(tag, attrs){ const n=document.createElementNS(SVGNS,tag); for(const k in attrs) n.setAttribute(k,attrs[k]); return n; }
function fmt(bp){ const n=Number(bp); return Number.isFinite(n)? n.toLocaleString("en-US") : String(bp); }
function num(v){ const n=Number(v); return Number.isFinite(n)?n:null; }

// ---- tooltip ----
const tip = $("#tip");
function showTip(e,html){ tip.innerHTML=html; tip.style.opacity=1; moveTip(e); }
function moveTip(e){ const p=14,w=tip.offsetWidth,h=tip.offsetHeight; let x=e.clientX+p,y=e.clientY+p;
  if(x+w>innerWidth)x=e.clientX-w-p; if(y+h>innerHeight)y=e.clientY-h-p; tip.style.left=x+"px"; tip.style.top=y+"px"; }
function hideTip(){ tip.style.opacity=0; }
function wire(node,html){ node.classList.add("hoverable");
  node.addEventListener("mouseenter",e=>showTip(e,html));
  node.addEventListener("mousemove",moveTip);
  node.addEventListener("mouseleave",hideTip); }

// ---- cohort picker ----
async function loadCohorts(){
  const sel = $("#cohort");
  try {
    const r = await fetch("/cohorts"); const d = await r.json();
    sel.innerHTML = "";
    const none = new Option("GENCODE only (auto-add enhancers)", ""); sel.add(none);
    let flagship = "";
    for(const c of d.cohorts){
      const o = new Option(`${c.name} · ${c.n_replicates} rep`, c.name);
      if(c.ontology_id === d.default) flagship = c.name;
      sel.add(o);
    }
    sel.value = ""; // default: let the server pick the flagship on enhancer questions
    sel.dataset.flagship = flagship;
  } catch(e){ /* offline cohort list; the empty option still works */ }
}

// ---- anchor: a chrN:pos(-pos) from the question or summary, for the variant marker ----
function parseAnchor(text){
  const m = /\b(chr[\w]+)\s*:\s*([\d,]+)(?:\s*-\s*([\d,]+))?/i.exec(text||"");
  if(!m) return null;
  const a = +m[2].replace(/,/g,""), b = m[3]? +m[3].replace(/,/g,"") : a;
  return { chrom:m[1], start:Math.min(a,b), end:Math.max(a,b), pos:Math.round((a+b)/2) };
}

// =================== render ===================
function renderResult(question, data){
  const main = $("#main"); main.innerHTML = "";

  const echo = document.createElement("p"); echo.className="qecho"; echo.textContent=question; main.appendChild(echo);

  const status = document.createElement("p"); status.className="status"; main.appendChild(status);
  if(data.note) status.innerHTML = `<span class="note">${escapeHtml(data.note)}</span>`;

  if(data.error){
    status.className = "status err"; status.textContent = data.error;
    if(data.code) main.appendChild(codePanel(data.code));
    return;
  }

  const summary = (data.summary||[]).join("\n");
  const records = data.records||[];
  const anchor = parseAnchor(question) || parseAnchor(summary);

  if(summary.trim()){
    const s = document.createElement("div"); s.className="summary mono"; s.textContent=summary; main.appendChild(s);
  }

  const plottable = records.filter(r => r.chrom!=null && num(r.start)!=null && num(r.end)!=null);
  if(plottable.length) main.appendChild(overviewBlock(plottable, anchor));

  const hasTx = records.some(r => (r.type||"")==="transcript");
  const hasExon = records.some(r => (r.type||"")==="exon" || (r.type||"")==="intron");
  if(anchor && (hasTx || hasExon)) main.appendChild(zoomBlock(records, anchor));

  if(data.code) main.appendChild(codePanel(data.code));
  if(records.length) main.appendChild(tableBlock(records));
  if(!plottable.length && !records.length && !summary.trim())
    { status.textContent = "(no result)"; }
}

function escapeHtml(s){ return String(s).replace(/[&<>"]/g, c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c])); }

function codePanel(code){
  const d = document.createElement("details"); d.className="code";
  const s = document.createElement("summary"); s.textContent="Generated pygenogrove code"; d.appendChild(s);
  const pre = document.createElement("pre"); pre.textContent=code; d.appendChild(pre);
  return d;
}

// ---- overview: gene models + enhancer lollipops + variant marker ----
function overviewBlock(recs, anchor){
  const genes = recs.filter(r => (r.type||"")==="gene");
  const enh   = recs.filter(r => (r.type||"")==="enhancer");
  const other = recs.filter(r => !["gene","enhancer","transcript","exon","intron"].includes(r.type||""));

  const block = section("Locus overview", enh.length ? "enhancer height = rE2G score" : "genomic coordinates");
  const W=1000, mL=62, mR=22;
  let lo=Infinity, hi=-Infinity;
  for(const r of recs){ lo=Math.min(lo,num(r.start)); hi=Math.max(hi,num(r.end)); }
  if(anchor){ lo=Math.min(lo,anchor.start); hi=Math.max(hi,anchor.end); }
  const span=Math.max(1,hi-lo), pad=span*0.03; const d0=lo-pad, d1=hi+pad;
  const x = bp => mL + (bp-d0)/(d1-d0)*(W-mL-mR);

  const laneY = { top:44 };
  const geneY0 = 50, geneRow = 20;
  const enhBase = enh.length ? 214 : 96, enhTop = 118;
  const bottom = enh.length ? enhBase : geneY0 + Math.max(1,genes.length+other.length)*geneRow + 8;
  const H = bottom + 20;
  const svg = el("svg", { viewBox:`0 0 ${W} ${H}`, role:"img" });

  // ruler
  ruler(svg, x, d0, d1, 26, W, mR, mL);

  // gene / generic feature models
  [...genes, ...other].forEach((g,i)=>{
    const y=geneY0+i*geneRow, gx0=x(num(g.start)), gx1=x(num(g.end));
    const isGene=(g.type||"")==="gene";
    const bar=el("rect",{x:gx0,y:y,width:Math.max(2,gx1-gx0),height:11,rx:2,
      fill:isGene?"var(--gene-fill)":"var(--feat)", stroke:isGene?"var(--gene)":"var(--feat)","stroke-width":1,
      opacity:isGene?1:.8});
    wire(bar, featTip(g)); svg.appendChild(bar);
    if(isGene && (g.strand==="+"||g.strand==="-")){
      const dir=g.strand==="+"?1:-1;
      for(let px=gx0+12; px<gx1-4; px+=26){
        svg.appendChild(el("path",{d:`M${px-3*dir},${y+3} L${px+3*dir},${y+5.5} L${px-3*dir},${y+8}`,
          fill:"none",stroke:"var(--gene)","stroke-width":1.1,"stroke-linecap":"round","stroke-linejoin":"round"}));
      }
    }
    const nm=el("text",{x:gx0,y:y-4,fill:"var(--ink2)","font-size":10.5,"font-weight":600}); nm.textContent=g.name||g.id||g.type||""; svg.appendChild(nm);
  });

  // enhancer score lane
  if(enh.length){
    [0,0.5,1].forEach(s=>{
      const y=enhBase - s*(enhBase-enhTop);
      svg.appendChild(el("line",{x1:mL,y1:y,x2:W-mR,y2:y,stroke:"var(--grid)","stroke-width":1,"stroke-dasharray":s?"2 3":"0"}));
      const t=el("text",{x:mL-8,y:y+3,fill:"var(--muted)","font-size":9.5,"text-anchor":"end"}); t.textContent=s.toFixed(1); svg.appendChild(t);
    });
    const axl=el("text",{x:mL-8,y:enhTop-8,fill:"var(--muted)","font-size":9.5,"text-anchor":"end"}); axl.textContent="score"; svg.appendChild(axl);
    for(const en of enh){
      const mid=x((num(en.start)+num(en.end))/2);
      const sc=num(en.score); const cls=en.class||"enhancer";
      const y=sc!=null ? enhBase - sc*(enhBase-enhTop) : enhTop;
      const col=CLS[cls]||"var(--feat)";
      svg.appendChild(el("line",{x1:mid,y1:enhBase,x2:mid,y2:y,stroke:col,"stroke-width":1.4,opacity:.55}));
      const dot=el("circle",{cx:mid,cy:y,r:4,fill:col,stroke:"var(--panel)","stroke-width":1.2});
      wire(dot, `<div class="t-h" style="color:${col}">${escapeHtml(cls)} enhancer${en.target?` → ${escapeHtml(en.target)}`:(en.gene?` → ${escapeHtml(en.gene)}`:"")}</div>`+
        (sc!=null?`<div class="t-r">rE2G score ${sc.toFixed(3)}${en.n?` · n=${en.n}`:""}</div>`:"")+
        `<div class="t-r">${escapeHtml(en.chrom||"")}:${fmt(en.start)}–${fmt(en.end)}</div>`);
      svg.appendChild(dot);
    }
  }

  // variant marker
  if(anchor){
    const vx=x(anchor.pos);
    svg.appendChild(el("line",{x1:vx,y1:laneY.top,x2:vx,y2:bottom,stroke:"var(--variant)","stroke-width":1.6}));
    svg.appendChild(el("path",{d:`M${vx},${laneY.top} l0,-8 l30,0 l0,13 l-30,0 z`,fill:"var(--variant)"}));
    const ft=el("text",{x:vx+4,y:laneY.top-2,fill:"#fff","font-size":9,"font-weight":700}); ft.textContent="VAR"; svg.appendChild(ft);
    const hit=el("rect",{x:vx-6,y:laneY.top,width:12,height:bottom-laneY.top,fill:"transparent"});
    wire(hit,`<div class="t-h" style="color:var(--variant)">Query position</div><div class="t-r">${escapeHtml(anchor.chrom)}:${fmt(anchor.pos)}</div>`);
    svg.appendChild(hit);
  }

  block.querySelector(".plot").appendChild(svg);
  block.appendChild(legend([
    anchor ? ['vline','var(--variant)','query position'] : null,
    genes.length ? ['bar','var(--gene)','gene model'] : null,
    other.length ? ['bar','var(--feat)','feature'] : null,
    enh.some(e=>(e.class)==="promoter") ? ['round','var(--promoter)','promoter'] : null,
    enh.some(e=>(e.class)==="genic") ? ['round','var(--genic)','genic'] : null,
    enh.some(e=>(e.class)==="intergenic") ? ['round','var(--intergenic)','intergenic'] : null,
  ]));
  return block;
}

// ---- zoom: exon vs intron across transcripts, at the anchor ----
function zoomBlock(recs, anchor){
  const txs = recs.filter(r => (r.type||"")==="transcript");
  const exons = recs.filter(r => (r.type||"")==="exon");
  const introns = recs.filter(r => (r.type||"")==="intron");

  // exon that sits over the anchor (its block defines the zoom window)
  const hitExons = exons.filter(e => num(e.start)<=anchor.pos && anchor.pos<=num(e.end));
  const exStart = hitExons.length ? Math.min(...hitExons.map(e=>num(e.start))) : anchor.pos-120;
  const exEnd   = hitExons.length ? Math.max(...hitExons.map(e=>num(e.end)))   : anchor.pos+120;
  const winPad = Math.max(140, (exEnd-exStart));
  const d0 = Math.min(exStart, anchor.start) - winPad, d1 = Math.max(exEnd, anchor.end) + winPad;

  // exon id looks like "exon:ENST....:20" — associate an exon to its transcript
  const txId = t => t.id || t.name;
  const exonTxOf = e => { const m=/exon:([^:]+):/.exec(e.id||""); return m?m[1]:null; };
  function callFor(t){
    const id = txId(t);
    const ex = hitExons.find(e => exonTxOf(e)===id) || hitExons.find(e => !exonTxOf(e));
    if(!ex) return { kind:"intron", coding:false };
    const cds = ex.cds;
    const coding = cds!=null && String(cds).toLowerCase()!=="none" && cds!=="";
    return { kind:"exon", coding, ex };
  }

  // rows: every transcript; if none, fall back to the exons themselves
  const rows = txs.length ? txs.map(t=>({t, call:callFor(t)}))
                          : exons.map(e=>({t:e, call:{kind:"exon", coding:(e.cds!=null&&String(e.cds).toLowerCase()!=="none"), ex:e}}));
  rows.sort((a,b)=> (a.call.kind===b.call.kind?0 : a.call.kind==="exon"?-1:1));

  const nExon = rows.filter(r=>r.call.kind==="exon").length;
  const block = section("At the query — exon vs. intron across transcripts",
                        `${nExon} exonic · ${rows.length-nExon} intronic`);
  const W=1000, mL=20, mR=170, rowH=26, y0=40;
  const H = y0 + rows.length*rowH + 20;
  const x = bp => mL + (bp-d0)/(d1-d0)*(W-mL-mR);
  const svg = el("svg",{viewBox:`0 0 ${W} ${H}`, role:"img"});

  // light grid every ~100bp
  const step = niceStep((d1-d0)/6);
  for(let t=Math.ceil(d0/step)*step; t<=d1; t+=step){
    const gx=x(t);
    svg.appendChild(el("line",{x1:gx,y1:24,x2:gx,y2:H-14,stroke:"var(--grid)","stroke-width":1}));
    const lab=el("text",{x:gx,y:16,fill:"var(--muted)","font-size":9.5,"text-anchor":"middle"}); lab.textContent=fmt(Math.round(t)); svg.appendChild(lab);
  }

  rows.forEach((row,i)=>{
    const {t,call}=row; const y=y0+i*rowH;
    const ix0=x(Math.max(d0,num(t.start ?? d0))), ix1=x(Math.min(d1,num(t.end ?? d1)));
    svg.appendChild(el("line",{x1:ix0,y1:y+7,x2:ix1,y2:y+7,stroke:"var(--axis)","stroke-width":2}));
    if(call.kind==="exon" && call.ex){
      svg.appendChild(el("rect",{x:x(num(call.ex.start)),y:y,width:Math.max(2,x(num(call.ex.end))-x(num(call.ex.start))),height:14,rx:2,
        fill: call.coding?"var(--exon)":"var(--exon-utr)"}));
    }
    const label = call.kind==="exon" ? (call.coding?"coding exon":"non-coding exon") : "intron";
    const col = call.kind==="exon" ? (call.coding?"var(--ink)":"var(--ink2)") : "var(--muted)";
    const id=el("text",{x:W-mR+12,y:y+6,fill:"var(--ink2)","font-size":10.5}); id.textContent=txId(t)||t.name||""; svg.appendChild(id);
    const cl=el("text",{x:W-mR+12,y:y+17,fill:col,"font-size":9.5}); cl.textContent=label; svg.appendChild(cl);
    const hit=el("rect",{x:mL,y:y-3,width:W-mR-mL,height:rowH-4,fill:"transparent"});
    wire(hit,`<div class="t-h">${escapeHtml(txId(t)||t.name||"")}</div>`+
      `<div class="t-r">${escapeHtml(t.chrom||anchor.chrom)}:${fmt(t.start)}–${fmt(t.end)}${t.strand?` (${t.strand})`:""}</div>`+
      `<div class="t-r">variant is <b>${label}</b> here</div>`);
    svg.appendChild(hit);
  });

  // variant marker
  const vx=x(anchor.pos);
  svg.appendChild(el("line",{x1:vx,y1:24,x2:vx,y2:H-14,stroke:"var(--variant)","stroke-width":1.6}));
  svg.appendChild(el("path",{d:`M${vx},24 l0,-9 l64,0 l0,13 l-64,0 z`,fill:"var(--variant)"}));
  const ft=el("text",{x:vx+5,y:22,fill:"#fff","font-size":9,"font-weight":700}); ft.textContent=fmt(anchor.pos); svg.appendChild(ft);

  block.querySelector(".plot").appendChild(svg);
  block.appendChild(legend([
    ['blk','var(--exon)','coding exon'],
    ['blk','var(--exon-utr)','non-coding exon'],
    ['bar','var(--axis)','intron'],
    ['vline','var(--variant)','variant position'],
  ]));
  return block;
}

function niceStep(raw){ const p=Math.pow(10,Math.floor(Math.log10(raw))); const f=raw/p;
  return (f>=5?5:f>=2?2:1)*p; }

function ruler(svg,x,d0,d1,yR,W,mR,mL){
  const step=niceStep((d1-d0)/6);
  for(let t=Math.ceil(d0/step)*step; t<=d1; t+=step){
    svg.appendChild(el("line",{x1:x(t),y1:yR,x2:x(t),y2:yR+5,stroke:"var(--axis)","stroke-width":1}));
    const lab=el("text",{x:x(t),y:yR-6,fill:"var(--muted)","font-size":10,"text-anchor":"middle"});
    lab.textContent = (d1-d0)>2e5 ? (t/1e6).toFixed(2)+" Mb" : fmt(Math.round(t)); svg.appendChild(lab);
  }
  svg.appendChild(el("line",{x1:mL,y1:yR,x2:W-mR,y2:yR,stroke:"var(--grid)","stroke-width":1}));
}

function section(title, sub){
  const s=document.createElement("section"); s.className="block";
  s.innerHTML=`<div class="cap"><h2>${escapeHtml(title)}</h2><p>${escapeHtml(sub||"")}</p></div><div class="plot"></div>`;
  return s;
}
function legend(items){
  const d=document.createElement("div"); d.className="legend";
  for(const it of items){ if(!it) continue; const [shape,color,label]=it;
    const cls = shape==="round"?"round":shape==="vline"?"vline":shape==="blk"?"blk":"";
    d.innerHTML += `<span><i class="${cls}" style="background:${color}${shape==="blk"?";height:9px;border-radius:2px":""}"></i> ${escapeHtml(label)}</span>`;
  }
  return d;
}

function featTip(g){
  const parts=[`<div class="t-h">${escapeHtml(g.name||g.id||g.type||"feature")}</div>`];
  if(g.id||g.biotype) parts.push(`<div class="t-r">${escapeHtml(g.id||"")}${g.biotype?" · "+escapeHtml(g.biotype):""}</div>`);
  parts.push(`<div class="t-r">${escapeHtml(g.chrom||"")}:${fmt(g.start)}–${fmt(g.end)}${g.strand?` (${g.strand})`:""}</div>`);
  return parts.join("");
}

// ---- records table (mirrors the CLI column order) ----
const COL_ORDER = ["chrom","start","end","strand","name","type","class","score","n","cohort","target","id","biotype"];
function tableBlock(records){
  const cols=[]; for(const k of COL_ORDER) if(records.some(r=>k in r)) cols.push(k);
  for(const r of records) for(const k in r) if(!cols.includes(k)) cols.push(k);
  const s=document.createElement("section"); s.className="block";
  s.innerHTML=`<div class="cap"><h2>Records</h2><p>${records.length} row${records.length===1?"":"s"}</p></div>`;
  const wrap=document.createElement("div"); wrap.className="tablewrap";
  const table=document.createElement("table");
  const head="<thead><tr>"+cols.map(c=>`<th>${escapeHtml(c)}</th>`).join("")+"</tr></thead>";
  const numCol=new Set(["start","end","score","n"]);
  const body="<tbody>"+records.map(r=>"<tr>"+cols.map(c=>{
    let v=r[c]; if(v==null) v="";
    if(c==="score" && num(v)!=null) v=num(v).toFixed(3);
    const sw = (c==="class"&&CLS[r[c]]) ? `<span class="swatch" style="background:${CLS[r[c]]}"></span>` : "";
    return `<td class="${numCol.has(c)?"mono":""}">${sw}${escapeHtml(v)}</td>`;
  }).join("")+"</tr>").join("")+"</tbody>";
  table.innerHTML=head+body; wrap.appendChild(table); s.appendChild(wrap);
  return s;
}

// ---- history ----
const history=[];
function addHistory(q){
  if(history[0]===q) return; history.unshift(q); if(history.length>30) history.pop();
  const ul=$("#hist"); ul.innerHTML="";
  for(const h of history){ const li=document.createElement("li"); li.textContent=h;
    li.title=h; li.onclick=()=>{ $("#q").value=h; submit(h); }; ul.appendChild(li); }
}

// ---- live progress steps ----
function addStep(msg){
  const ul=$("#steps"); if(!ul) return;
  const prev=ul.lastElementChild;
  if(prev){ prev.classList.remove("active"); prev.classList.add("done"); const t=prev.querySelector(".tick"); if(t) t.textContent=""; }
  const li=document.createElement("li"); li.className="active";
  li.innerHTML=`<span class="spin"></span><span class="lbl">${escapeHtml(msg)}</span><span class="tick"></span>`;
  ul.appendChild(li);
}
function tickStep(msg){ const li=$("#steps")?.lastElementChild; const t=li?.querySelector(".tick"); if(t) t.textContent=msg; }

// ---- submit: stream ndjson events, narrate steps, render on done ----
let busy=false;
async function submit(question){
  question=(question||"").trim(); if(!question||busy) return;
  busy=true; const btn=$(".ask"); btn.disabled=true; btn.textContent="…";
  const main=$("#main");
  main.innerHTML=`<p class="qecho">${escapeHtml(question)}</p><div class="progress"><ul class="steps" id="steps"></ul></div>`;
  addStep("Sending the question…");
  try{
    const r=await fetch("/ask",{method:"POST",headers:{"Content-Type":"application/json"},
      body:JSON.stringify({question, cohort:$("#cohort").value})});
    if(!r.ok){ throw new Error((await r.text())||("HTTP "+r.status)); }
    const reader=r.body.getReader(), dec=new TextDecoder(); let buf="", final=null;
    for(;;){
      const {value,done}=await reader.read(); if(done) break;
      buf+=dec.decode(value,{stream:true}); let i;
      while((i=buf.indexOf("\n"))>=0){
        const line=buf.slice(0,i).trim(); buf=buf.slice(i+1);
        if(!line) continue;
        const ev=JSON.parse(line);
        if(ev.t==="step") addStep(ev.msg);
        else if(ev.t==="tick") tickStep(ev.msg);
        else if(ev.t==="done") final=ev.result;
        else if(ev.t==="error") final={error:ev.msg};
      }
    }
    renderResult(question, final||{error:"(no response from server)"}); addHistory(question);
  }catch(e){
    main.innerHTML=`<p class="qecho">${escapeHtml(question)}</p><p class="status err">Request failed: ${escapeHtml(String(e))}</p>`;
  }finally{ busy=false; btn.disabled=false; btn.textContent="Ask"; }
}

$("#askform").addEventListener("submit", e=>{ e.preventDefault(); submit($("#q").value); });
document.querySelectorAll(".placeholder code").forEach(c=>{
  c.onclick=()=>{ $("#q").value=c.textContent; submit(c.textContent); };
});
loadCohorts();
$("#q").focus();
</script>
</body>
</html>
"""
