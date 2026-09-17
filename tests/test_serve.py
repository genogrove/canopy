# SPDX-License-Identifier: GPL-3.0-or-later
"""Backend checks for ``canopy serve`` — cohort fallback + request routing.

The record→SVG rendering is client-side JS (exercised in a browser); here we guard the
Python glue: an enhancer question with no picker choice falls back to the flagship cohort,
a structural question stays on plain GENCODE, and /ask returns the record contract while a
bad request is rejected — all without an API key or a grove (the pipeline is stubbed)."""

import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from genogrove_canopy import serve


def test_serve_cohort_precedence():
    E = ["enhancers"]
    # UI picker overrides the model's declaration
    cohorts, note = serve._serve_cohorts("LNCaP", "K562", E)
    assert cohorts and note is None
    # else the model's declared cohort
    cohorts, note = serve._serve_cohorts("", "K562", E)
    assert cohorts and note is None
    # neither -> the default cohort (flagged) — for an enhancer question only
    cohorts, note = serve._serve_cohorts("", "", E)
    assert cohorts and note == "default"
    # an SV question with no cohort gets no default: SVs are per tumour cohort
    cohorts, note = serve._serve_cohorts("", "", ["sv"])
    assert cohorts == {} and "tumour cohort" in note
    # a declared tissue with no catalog match -> none, with a note (never silently substituted)
    cohorts, note = serve._serve_cohorts("", "NoSuchTissue12345", E)
    assert cohorts == {} and note


def _fake_pipeline(question, cohort, model, emit):
    emit("step", "Opening the grove")
    emit("tick", "fetched 2/2 replicate tracks")
    emit("step", "Running the query over the grove")
    return {"summary": [f"q: {question}"], "records": [{"chrom": "chr7", "start": 1, "end": 9}],
            "code": "pass", "note": "note"}


@pytest.fixture
def server(monkeypatch):
    monkeypatch.setattr(serve, "_pipeline", _fake_pipeline)
    srv = ThreadingHTTPServer(("127.0.0.1", 0), serve._Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


def _events(url, payload):
    """POST /ask and collect the streamed ndjson events."""
    req = urllib.request.Request(url + "/ask", data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    body = urllib.request.urlopen(req, timeout=5).read().decode()
    return [json.loads(line) for line in body.splitlines() if line.strip()]


def test_ask_streams_steps_then_result(server):
    events = _events(server, {"question": "x", "cohort": ""})
    kinds = [e["t"] for e in events]
    assert "step" in kinds and "tick" in kinds        # progress was narrated
    assert kinds[-1] == "done"                         # result comes last
    result = events[-1]["result"]
    assert result["records"] and result["code"] == "pass"


def test_empty_question_rejected(server):
    req = urllib.request.Request(server + "/ask", data=json.dumps({"question": "", "cohort": ""}).encode(),
                                 headers={"Content-Type": "application/json"})
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(req, timeout=5)
    assert e.value.code == 400


@pytest.mark.parametrize("body,headers", [
    (b"{not json", {"Content-Type": "application/json"}),   # unparseable body
    (b'{"question":"x"}', {"Content-Length": "abc"}),        # non-integer Content-Length
    (b'{"question": 123}', {}),                              # question isn't a string
])
def test_malformed_request_gets_400(server, body, headers):
    """Each of these used to raise before send_response, so the client got no reply at all
    (RemoteDisconnected) and a traceback hit stderr — not the 400 this endpoint promises."""
    req = urllib.request.Request(server + "/ask", data=body, headers=headers)
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(req, timeout=5)
    assert e.value.code == 400


def test_cohort_default_flag_is_not_display_text(monkeypatch):
    """A cohort that resolves but yields no enhancer links must not surface the internal
    "default" sentinel as the user-facing note (the CLI filters it; serve must too)."""
    monkeypatch.setattr(serve, "_resolve_cohorts",
                        lambda specs: {"LNCaP clone FGC": {"re2g": ["EFO:0005726"], "pcawg": []}})
    monkeypatch.setattr(serve.resources, "_all_grove_gg",
                        lambda n: type("P", (), {"exists": lambda s: True})())
    monkeypatch.setattr(serve, "_grove", lambda m: type("G", (), {
        "system_prompt": "", "preamble": "", "gg": "", "model": m,
        "worker": type("W", (), {"submit": lambda s, c: type("R", (), {
            "returncode": 0, "timed_out": False, "truncated": False, "stdout": "ok: 1", "stderr": ""})()})()})())
    monkeypatch.setattr(serve.llm, "generate_query",
                        lambda q, sp, model=None: ("", ["enhancers"], "pass"))
    from genogrove_canopy.layers import enhancers
    monkeypatch.setattr(enhancers, "ensure_index", lambda cohort: False)

    result = serve._pipeline("enhancers of MYC", "", "m", lambda k, m: None)
    assert result["note"] != "default" and "no enhancers attached" in result["note"]   # the gap is said, never the sentinel


def test_page_served(server):
    body = urllib.request.urlopen(server + "/", timeout=5).read()
    assert body.lstrip().startswith(b"<!doctype")
