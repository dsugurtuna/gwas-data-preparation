"""Strand classification and merge planning, with PLINK replaced by a fake."""

from pathlib import Path

import pytest

from gwas_prep.assembler import GenotypeAssembler, compare_alleles


@pytest.mark.parametrize(
    ("ref", "other", "expected"),
    [
        (("A", "G"), ("A", "G"), "same"),
        (("A", "G"), ("G", "A"), "same"),
        (("A", "G"), ("T", "C"), "flipped"),
        (("A", "G"), ("C", "T"), "flipped"),
        (("A", "T"), ("A", "T"), "ambiguous"),
        (("C", "G"), ("G", "C"), "ambiguous"),
        (("A", "G"), ("A", "C"), "incompatible"),
        (("0", "G"), ("A", "G"), "same"),
        (("I", "D"), ("D", "I"), "same"),
    ],
)
def test_compare_alleles(
    ref: tuple[str, str], other: tuple[str, str], expected: str
) -> None:
    assert compare_alleles(ref, other) == expected


def _bim(path: Path, rows: list[tuple[str, str, str]]) -> None:
    path.write_text(
        "".join(
            f"1 {vid} 0 {i + 1}00 {a1} {a2}\n" for i, (vid, a1, a2) in enumerate(rows)
        )
    )


class FakePlink:
    """Records commands; writes a .bim/.fam for every --make-bed output."""

    def __init__(self, fail_on: str | None = None) -> None:
        self.calls: list[list[str]] = []
        self.fail_on = fail_on

    def __call__(self, argv: list[str]) -> int:
        self.calls.append(argv)
        if self.fail_on and self.fail_on in argv:
            out = argv[argv.index("--out") + 1]
            Path(f"{out}-merge.missnp").write_text("rs9\nrs10\n")
            return 3
        if "--make-bed" in argv:
            out = argv[argv.index("--out") + 1]
            Path(f"{out}.bim").write_text("1 rs1 0 1 A G\n1 rs2 0 2 C T\n")
            Path(f"{out}.fam").write_text(
                "F1 S1 0 0 1 -9\nF2 S2 0 0 2 -9\nF3 S3 0 0 1 -9\n"
            )
        return 0


@pytest.fixture()
def batches(tmp_path: Path) -> list[Path]:
    b1, b2, b3 = tmp_path / "b1", tmp_path / "b2", tmp_path / "b3"
    _bim(
        Path(f"{b1}.bim"),
        [("rs1", "A", "G"), ("rs2", "C", "T"), ("rs3", "A", "T"), ("rs4", "A", "G")],
    )
    _bim(
        Path(f"{b2}.bim"),
        [("rs1", "T", "C"), ("rs2", "C", "T"), ("rs3", "T", "A"), ("rs4", "A", "C")],
    )
    _bim(Path(f"{b3}.bim"), [("rs1", "A", "G"), ("rs2", "C", "T")])
    return [b1, b2, b3]


def test_check_strand(batches: list[Path]) -> None:
    check = GenotypeAssembler.check_strand(f"{batches[0]}.bim", f"{batches[1]}.bim")
    assert check.flipped == ["rs1"]
    assert check.ambiguous == ["rs3"]
    assert check.incompatible == ["rs4"]


def test_detect_strand_conflicts_space_delimited(batches: list[Path]) -> None:
    conflicts = GenotypeAssembler().detect_strand_conflicts(
        f"{batches[0]}.bim", f"{batches[1]}.bim"
    )
    assert conflicts == ["rs1", "rs3"]


def test_merge_flips_only_where_needed(batches: list[Path], tmp_path: Path) -> None:
    fake = FakePlink()
    ga = GenotypeAssembler(work_dir=tmp_path / "work", runner=fake)
    result = ga.merge_batches(batches, tmp_path / "merged")
    assert result.success
    flips = [c for c in fake.calls if "--flip" in c]
    assert len(flips) == 1 and flips[0][2] == str(batches[1])  # only batch 2
    assert (tmp_path / "work" / "batch1.flip").read_text() == "rs1\n"
    assert result.strand_flips == 1
    assert result.ambiguous_variants == 1
    # rs4 is incompatible and is excluded from every batch.
    assert result.excluded_variants == 1
    assert sum("--exclude" in c for c in fake.calls) == 3
    assert result.total_samples == 3
    assert result.total_variants == 2


def test_merge_can_exclude_ambiguous(batches: list[Path], tmp_path: Path) -> None:
    ga = GenotypeAssembler(work_dir=tmp_path / "work", runner=FakePlink())
    result = ga.merge_batches(batches, tmp_path / "merged", exclude_ambiguous=True)
    assert result.excluded_variants == 2
    assert (tmp_path / "work" / "exclude.txt").read_text() == "rs3\nrs4\n"


def test_merge_failure_reports_missnp(batches: list[Path], tmp_path: Path) -> None:
    ga = GenotypeAssembler(
        work_dir=tmp_path / "work", runner=FakePlink(fail_on="--merge-list")
    )
    result = ga.merge_batches(batches, tmp_path / "merged")
    assert not result.success
    assert result.error is not None and "2 variants listed" in result.error


def test_count_with_dotted_prefix(tmp_path: Path) -> None:
    Path(f"{tmp_path}/cohort.chr1.fam").write_text("F S 0 0 1 -9\n")
    assert GenotypeAssembler.count_samples(tmp_path / "cohort.chr1") == 1
