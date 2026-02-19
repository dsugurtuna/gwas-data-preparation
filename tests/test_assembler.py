"""Tests for GenotypeAssembler."""

from gwas_prep.assembler import GenotypeAssembler, AssemblyResult


class TestGenotypeAssembler:
    def test_count_samples(self, tmp_path):
        fam = tmp_path / "test.fam"
        fam.write_text("FAM1 IID1 0 0 1 1\nFAM2 IID2 0 0 2 1\nFAM3 IID3 0 0 1 2\n")
        ga = GenotypeAssembler()
        assert ga.count_samples(tmp_path / "test") == 3

    def test_count_variants(self, tmp_path):
        bim = tmp_path / "test.bim"
        bim.write_text("1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tC\tT\n")
        ga = GenotypeAssembler()
        assert ga.count_variants(tmp_path / "test") == 2

    def test_count_missing_file(self, tmp_path):
        ga = GenotypeAssembler()
        assert ga.count_samples(tmp_path / "missing") == 0
        assert ga.count_variants(tmp_path / "missing") == 0

    def test_detect_strand_conflicts(self, tmp_path):
        bim_a = tmp_path / "a.bim"
        bim_a.write_text("1\trs1\t0\t100\tA\tT\n1\trs2\t0\t200\tC\tG\n1\trs3\t0\t300\tA\tG\n")
        bim_b = tmp_path / "b.bim"
        bim_b.write_text("1\trs1\t0\t100\tT\tA\n1\trs2\t0\t200\tG\tC\n1\trs3\t0\t300\tA\tG\n")
        ga = GenotypeAssembler()
        conflicts = ga.detect_strand_conflicts(bim_a, bim_b)
        assert "rs1" in conflicts
        assert "rs2" in conflicts
        assert "rs3" not in conflicts

    def test_merge_too_few_batches(self):
        ga = GenotypeAssembler()
        result = ga.merge_batches(["single"], "output")
        assert not result.success
        assert "two batches" in result.error


class TestAssemblyResult:
    def test_default_success(self):
        r = AssemblyResult()
        assert r.success is True
        assert r.batches_merged == 0
