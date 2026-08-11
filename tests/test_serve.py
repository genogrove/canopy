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

from canopy import serve


def test_serve_cohort_precedence():
    # UI picker overrides the model's declaration
    cohorts, note = serve._serve_cohorts("LNCaP", "K562")
    assert cohorts and note is None
    # else the model's declared cohort
    cohorts, note = serve._serve_cohorts("", "K562")
    assert cohorts and note is None
    # neither -> the default cohort (flagged)
    cohorts, note = serve._serve_cohorts("", "")
    assert cohorts and note == "default"
    # a declared tissue with no catalog match -> none, with a note (never silently substituted)
    cohorts, note = serve._serve_cohorts("", "NoSuchTissue12345")
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


def test_page_served(server):
    body = urllib.request.urlopen(server + "/", timeout=5).read()
    assert body.lstrip().startswith(b"<!doctype")
