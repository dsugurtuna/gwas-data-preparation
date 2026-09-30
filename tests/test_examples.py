"""The synthetic examples contain known, planted failures; check they are found."""

import runpy
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_qc_demo_finds_planted_failures(capsys: pytest.CaptureFixture[str]) -> None:
    runpy.run_path(str(EXAMPLES / "qc_demo.py"), run_name="__main__")
    out = capsys.readouterr().out
    assert "variant call_rate       ['rs100007']" in out
    assert "variant hwe             ['rs100033']" in out
    assert "sample heterozygosity  ['SYN027']" in out
    assert "sample relatedness     ['SYN004']" in out
    assert "variants kept 56/60, samples kept 36/40" in out
    assert "flipped ['rs1'], ambiguous ['rs3'], incompatible ['rs4']" in out
