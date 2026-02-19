"""Tests for QualityController."""

from gwas_prep.qc import QualityController, QCReport


class TestQualityController:
    def test_check_variant_call_rates(self, tmp_path):
        lmiss = tmp_path / "data.lmiss"
        lmiss.write_text(
            "CHR SNP N_MISS N_GENO F_MISS\n"
            "1 rs1 2 100 0.02\n"
            "1 rs2 5 100 0.05\n"
            "1 rs3 1 100 0.01\n"
        )
        qc = QualityController(call_rate_variant=0.98)
        failed = qc.check_variant_call_rates(lmiss)
        assert "rs2" in failed
        assert "rs1" not in failed
        assert "rs3" not in failed

    def test_check_sample_call_rates(self, tmp_path):
        imiss = tmp_path / "data.imiss"
        imiss.write_text(
            "FID IID MISS_PHENO N_MISS N_GENO F_MISS\n"
            "F1 S1 N 10 1000 0.01\n"
            "F2 S2 N 50 1000 0.05\n"
        )
        qc = QualityController(call_rate_sample=0.98)
        failed = qc.check_sample_call_rates(imiss)
        assert "S2" in failed
        assert "S1" not in failed

    def test_check_maf(self, tmp_path):
        frq = tmp_path / "data.frq"
        frq.write_text(
            "CHR SNP A1 A2 MAF NCHROBS\n"
            "1 rs1 A G 0.15 200\n"
            "1 rs2 C T 0.005 200\n"
        )
        qc = QualityController(maf_threshold=0.01)
        failed = qc.check_maf(frq)
        assert "rs2" in failed
        assert "rs1" not in failed


class TestQCReport:
    def test_generate_report(self):
        qc = QualityController()
        report = qc.generate_report(
            initial_samples=1000,
            initial_variants=50000,
            failed_variants={
                "call_rate": {"rs1", "rs2"},
                "maf": {"rs3"},
                "hwe": {"rs4", "rs5"},
            },
            failed_samples={
                "call_rate": {"S1"},
                "heterozygosity": {"S2", "S3"},
            },
        )
        assert report.final_variants == 50000 - 5
        assert report.final_samples == 1000 - 3
        assert report.removed_low_call_rate_variants == 2
        assert report.removed_low_maf_variants == 1
        assert report.removed_hwe_variants == 2

    def test_pass_rates(self):
        r = QCReport(initial_samples=100, initial_variants=1000, final_samples=90, final_variants=950)
        assert r.sample_pass_rate == 0.9
        assert r.variant_pass_rate == 0.95

    def test_zero_division(self):
        r = QCReport()
        assert r.sample_pass_rate == 0.0
        assert r.variant_pass_rate == 0.0
