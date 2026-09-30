"""FormatConverter commands, with PLINK replaced by a fake."""

import gzip
import subprocess
from pathlib import Path

from gwas_prep.converter import FormatConverter

VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS_1\tS_2\n"
    "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n"
)


class FakePlink:
    def __init__(self) -> None:
        self.argv: list[str] = []

    def __call__(self, argv: list[str]) -> subprocess.CompletedProcess[str]:
        self.argv = argv
        out = argv[argv.index("--out") + 1]
        if "--make-bed" in argv:
            Path(f"{out}.fam").write_text("S_1 S_1 0 0 0 -9\nS_2 S_2 0 0 0 -9\n")
            Path(f"{out}.bim").write_text("1\trs1\t0\t100\tG\tA\n")
        elif "bgz" in argv:
            with gzip.open(f"{out}.vcf.gz", "wt") as fh:
                fh.write(VCF)
        else:
            Path(f"{out}.vcf").write_text(VCF)
        return subprocess.CompletedProcess(argv, 0, "", "")


def test_vcf_to_plink_keeps_whole_sample_ids(tmp_path: Path) -> None:
    fake = FakePlink()
    result = FormatConverter(runner=fake).vcf_to_plink("in.vcf", tmp_path / "out")
    assert "--double-id" in fake.argv
    assert "--keep-allele-order" in fake.argv
    assert (result.sample_count, result.variant_count) == (2, 1)


def test_plink_to_vcf_uses_iid_and_keeps_ref(tmp_path: Path) -> None:
    fake = FakePlink()
    result = FormatConverter(runner=fake).plink_to_vcf("in", tmp_path / "out.vcf")
    assert fake.argv[3:5] == ["--recode", "vcf-iid"]
    assert "--keep-allele-order" in fake.argv
    assert fake.argv[-1] == str(tmp_path / "out")  # PLINK appends .vcf itself
    assert result.output_path == str(tmp_path / "out.vcf")
    assert (result.sample_count, result.variant_count) == (2, 1)


def test_plink_to_vcf_gz(tmp_path: Path) -> None:
    fake = FakePlink()
    result = FormatConverter(runner=fake).plink_to_vcf("in", tmp_path / "out.vcf.gz")
    assert "bgz" in fake.argv
    assert result.output_path.endswith("out.vcf.gz")
    assert result.sample_count == 2


def test_failure_is_reported(tmp_path: Path) -> None:
    def failing(argv: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(argv, 2, "", "Error: No '_' in sample ID.")

    result = FormatConverter(runner=failing).vcf_to_plink("in.vcf", tmp_path / "o")
    assert not result.success
    assert result.error == "Error: No '_' in sample ID."
