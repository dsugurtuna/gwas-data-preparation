"""Tests for QualityController."""

from gwas_prep.qc import QCReport, QualityController


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
            "CHR SNP A1 A2 MAF NCHROBS\n1 rs1 A G 0.15 200\n1 rs2 C T 0.005 200\n"
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
        r = QCReport(
            initial_samples=100,
            initial_variants=1000,
            final_samples=90,
            final_variants=950,
        )
        assert r.sample_pass_rate == 0.9
        assert r.variant_pass_rate == 0.95

    def test_zero_division(self):
        r = QCReport()
        assert r.sample_pass_rate == 0.0
        assert r.variant_pass_rate == 0.0


def _write(path, text):
    path.write_text(text)
    return path


def test_check_hwe_uses_requested_test_rows(tmp_path):
    hwe = _write(
        tmp_path / "d.hwe",
        " CHR  SNP  TEST  A1 A2  GENO  O(HET)  E(HET)  P\n"
        "   1  rs1  ALL    A  G  1/2/97  0.02  0.03  1e-08\n"
        "   1  rs1  UNAFF  A  G  1/2/47  0.04  0.04  0.5\n"
        "   1  rs2  ALL    C  T  30/40/30  0.4  0.42  0.2\n",
    )
    qc = QualityController(hwe_p=1e-6)
    assert qc.check_hwe(hwe) == {"rs1"}
    assert qc.check_hwe(hwe, test="UNAFF") == set()


def test_check_heterozygosity_outliers(tmp_path):
    rows = "".join(f"F{i} S{i} 700 700 1000 0.0{i % 3}\n" for i in range(30))
    het = _write(
        tmp_path / "d.het",
        "FID IID O(HOM) E(HOM) N(NM) F\n" + rows + "FX SX 600 700 1000 -0.40\n",
    )
    assert QualityController(het_sd=3.0).check_heterozygosity(het) == {"SX"}


def test_check_sex(tmp_path):
    sexcheck = _write(
        tmp_path / "d.sexcheck",
        "FID IID PEDSEX SNPSEX STATUS F\nF1 S1 1 1 OK 0.99\nF2 S2 2 1 PROBLEM 0.97\n",
    )
    assert QualityController.check_sex(sexcheck) == {"S2"}


def test_check_relatedness_removes_fewest(tmp_path):
    header = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n"
    rows = [
        ("A", "B", 0.5),
        ("A", "C", 0.5),
        ("A", "D", 0.25),
        ("E", "F", 0.1),
    ]
    body = "".join(f"x {a} x {b} UN NA 0 1 0 {p} -1 0.9 1 NA\n" for a, b, p in rows)
    genome = _write(tmp_path / "d.genome", header + body)
    # Removing A alone breaks every pair above 0.2; E/F is below threshold.
    assert QualityController(pi_hat=0.2).check_relatedness(genome) == {"A"}


def test_maf_na_counts_as_failed(tmp_path):
    frq = _write(tmp_path / "d.frq", "CHR SNP A1 A2 MAF NCHROBS\n1 rs1 0 G NA 200\n")
    assert QualityController().check_maf(frq) == {"rs1"}


def test_report_counts_union_once():
    qc = QualityController()
    report = qc.generate_report(
        100,
        1000,
        {"call_rate": {"rs1", "rs2"}, "maf": {"rs2", "rs3"}},
        {"sex": {"S1"}, "relatedness": {"S1", "S2"}},
    )
    assert report.final_variants == 997
    assert report.final_samples == 98
