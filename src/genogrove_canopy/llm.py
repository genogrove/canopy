# SPDX-License-Identifier: GPL-3.0-or-later
"""Code generation via the Anthropic API.

Given a natural-language question, ask Claude to emit Python that drives
``pygenogrove`` to compute the answer. The model is given the ``pygenogrove`` API
surface, the grove model, and the curated resource context (see
:mod:`genogrove_canopy.resources`) via the system prompt in ``prompts/system.md``.

The generated code is untrusted: it must only ever run through :mod:`genogrove_canopy.sandbox`,
never ``exec``'d here.
"""

from __future__ import annotations

import re
from pathlib import Path

DEFAULT_MODEL = "claude-opus-4-8"

_SYSTEM_MD = Path(__file__).with_name("prompts") / "system.md"

#: Ways a model spells "no biosample" instead of omitting the line as system.md asks. Normalised
#: to "" so the host doesn't hand a placeholder to the ENCODE catalog. This is deliberately
#: permissive: canopy is meant to be model-agnostic, and treating one vendor's phrasing as the
#: only valid one would score every other model down for a formatting habit, not a capability.
_NO_COHORT = frozenset((
    "", "none", "n/a", "na", "null", "nil", "-", "--",
    "unspecified", "not specified", "not applicable",
))
# The system prompt's runtime-injected datasets block (the TODO placeholder).
_RESOURCES_HEADING = "## Available resources"


def build_system_prompt(resources_block: str) -> str:
    """The codegen system prompt: ``system.md`` with the resources block injected.

    ``resources_block`` replaces the placeholder under "Available resources" — it
    names the variables holding each dataset path and what they are. The generated
    code emits canonical records (see the Rules); the host renders the user's
    ``--format`` choice, so the format is not part of this prompt.
    """
    text = _SYSTEM_MD.read_text(encoding="utf-8")
    head, _, _tail = text.partition(_RESOURCES_HEADING)
    return f"{head}{_RESOURCES_HEADING}\n\n{resources_block.strip()}\n"


def generate_query(question: str, system_prompt: str, *, model: str = DEFAULT_MODEL, on_usage=None):
    """Translate ``question`` into ``pygenogrove`` Python via Claude.

    Returns ``(cohort, targets, code)``: ``cohort`` is the biosample/cell-line(s) the model read
    from the question (``""`` if none named; ``;``-separated if several; the host resolves it
    against the ENCODE catalog and attaches those cohorts into the grove before running
    ``code``), ``targets`` is the declared enhancer targets — ``{"gene": ...}`` /
    ``{"region": ...}`` — whose presence marks the question as regulatory (empty otherwise).
    ``code`` is the generated Python. The caller runs it through the
    sandbox; nothing is executed here. Raises ``RuntimeError`` if the model declines.

    ``on_usage``, if given, is called with the response's ``usage`` object before parsing — a
    side channel for cost accounting (a benchmark harness, not part of the package) that keeps the return shape stable for
    ``cli``/``serve``, which don't care. A callback rather than module state because ``serve``
    answers questions on multiple threads.
    """
    import anthropic  # lazy: keeps the module importable without the SDK/key

    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the environment
    response = client.messages.create(
        model=model,
        max_tokens=16000,  # room for adaptive thinking + a small program; non-streaming-safe
        thinking={"type": "adaptive"},
        system=system_prompt,
        messages=[{"role": "user", "content": question}],
    )
    if on_usage is not None:
        on_usage(response.usage)
    if response.stop_reason == "refusal":
        raise RuntimeError("the model declined to answer this question")
    text = "".join(b.text for b in response.content if b.type == "text")
    return parse_targets_and_code(text)


def parse_targets_and_code(text: str):
    """Split the model's reply into ``(cohort, targets, code)``.

    ``cohort`` comes from an optional ``COHORT: <name>`` line, ``targets`` from an optional
    ``TARGETS: [ ... ]`` JSON line (both outside the code fence; absent/malformed → ``""`` / ``[]``).
    ``code`` is the fenced program (or the text as-is)."""
    import json

    cohort = ""
    # Horizontal whitespace only: `\s` matches newlines, so `COHORT:\s*(.+?)` would step over a
    # bare `COHORT:` line (the correct declaration for a question naming no biosample) and capture
    # the *following* line — reading the cohort as "TARGETS:" and sending that to the catalog.
    mc = re.search(r"^[ \t]*COHORT:[ \t]*(.+?)[ \t]*$", text, re.MULTILINE)
    if mc:
        cohort = mc.group(1).strip().strip("\"'").strip("()").strip()
        if cohort.lower() in _NO_COHORT:
            cohort = ""

    targets = []
    mt = re.search(r"^\s*TARGETS:\s*(\[.*?\])\s*$", text, re.MULTILINE | re.DOTALL)
    if mt:
        try:
            parsed = json.loads(mt.group(1))
            if isinstance(parsed, list):
                targets = [t for t in parsed if isinstance(t, dict)]
        except ValueError:
            pass

    code = _strip_code_fence(text)
    # Defensive: strip any COHORT:/TARGETS: declaration lines that leaked into the program body
    # (e.g. the model fenced them with the code). `COHORT: MCF-7` is a valid Python annotation
    # that would NameError at runtime, so never let it reach the sandbox.
    code = "\n".join(ln for ln in code.splitlines()
                     if not re.match(r"\s*(COHORT|TARGETS)\s*:", ln)).strip() + "\n"
    return cohort, targets, code


def _strip_code_fence(text: str) -> str:
    """Return the Python program from the model's reply.

    Prefer an explicit ```python fence (so a bare ``` block holding the COHORT/TARGETS
    declarations is never mistaken for the program); else the last ``` block; else the text
    as-is. system.md asks for one ```python fence, but be robust to variation.
    """
    m = re.search(r"```python\n(.*?)```", text, re.DOTALL)
    if not m:
        blocks = re.findall(r"```(?:\w+)?\n(.*?)```", text, re.DOTALL)
        return (blocks[-1] if blocks else text).strip() + "\n"
    return m.group(1).strip() + "\n"
