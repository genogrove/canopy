# SPDX-License-Identifier: GPL-3.0-or-later
"""Checks on the codegen contract that need **no bindings**, so they run in the default CI job.

The counterpart, ``test_api_surface.py``, asserts the advertised methods exist on the real
``GroveView`` and therefore only binds in the ``api-surface`` workflow. Everything here is pure
Python over a tuple and a rendered string — keep it that way, or these stop running on most PRs.
"""

from __future__ import annotations

from genogrove_canopy.cli import QUERY_SURFACE, resources_block


def test_advertised_surface_names_no_known_mutator() -> None:
    """None of the advertised names is a known grove **writer**.

    A lint on the constant, not a check against the build: it compares ``QUERY_SURFACE`` against a
    hardcoded list of mutating method names. It cannot detect a future ``pygenogrove`` that moves a
    writer onto ``GroveView`` — that would leave ``QUERY_SURFACE`` untouched and this test passing.
    Proving a method does not mutate means calling it, and calling writers against a real grove to
    see whether they throw is a worse test than none. What this does catch is someone extending the
    "Query-only" list with a name that plainly is not.
    """
    mutating = {"insert", "insert_bulk", "insert_sorted", "add_edge", "remove_edge", "remove_key",
                "remove_edges_from", "remove_edges_if", "serialize", "compact", "clear_graph"}
    overlap = sorted(set(QUERY_SURFACE) & mutating)
    assert not overlap, f"QUERY_SURFACE advertises mutating method(s) {overlap} as query-only"


def test_every_advertised_method_reaches_the_model() -> None:
    """The advertised names must appear in the **rendered** resources block.

    Checks the contract rather than the source text: an earlier version asserted the token
    ``QUERY_SURFACE`` appeared somewhere in ``cli.py`` after its definition, which a comment
    mention would have satisfied. Rendering the block is what proves the names actually reach the
    model.
    """
    assert QUERY_SURFACE, "QUERY_SURFACE is empty — the prompt would advertise no query surface"

    block = resources_block("GENCODE_HUMAN", "test description", "- test layer")
    for method in QUERY_SURFACE:
        assert f"`{method}`" in block, f"{method} is in QUERY_SURFACE but not in the rendered block"
    assert "Read-only" in block


def test_prompt_names_only_allowlisted_modules() -> None:
    """Every module the system prompt imports in a code fence, or names in backticks as a thing
    to use, must be importable in the sandbox. The prompt once told the model to parse own files
    "with `csv`" — a module the allowlist excludes by name — so the first own-data turn burned
    on ``ImportError``. The prompt is the contract; the allowlist is the enforcement; they must
    agree."""
    import importlib.util
    import re
    from pathlib import Path

    from genogrove_canopy import llm, sandbox

    md = Path(llm.__file__).with_name("prompts").joinpath("system.md").read_text()
    imported = {m.group(1).split(".")[0]
                for m in re.finditer(r"^\s*(?:import|from)\s+([\w.]+)", md, re.M)}
    # a backticked word that resolves to an importable module is a module the prompt recommends
    named = {w for w in set(re.findall(r"`(\w+)`", md))
             if w not in ("GROVE",) and importlib.util.find_spec(w) is not None}
    bad = sorted((imported | named) - set(sandbox.ALLOWED_IMPORTS))
    assert not bad, f"system.md names modules the sandbox blocks: {bad}"


def test_resources_block_interpolates_its_arguments() -> None:
    """Guards the inverse: a block that ignored its inputs would make the test above vacuous."""
    block = resources_block("MY_VAR", "MY_DESCRIPTION", "MY_LAYERS")
    assert "MY_VAR" in block and "MY_DESCRIPTION" in block and "MY_LAYERS" in block
