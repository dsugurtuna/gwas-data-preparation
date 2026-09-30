"""PLINK binary <-> VCF conversion with PLINK 1.9.

Two PLINK defaults are overridden on purpose:

* On import, PLINK 1.9 splits a VCF sample ID at its single underscore into
  FID and IID (``SAMPLE_001`` becomes FID ``SAMPLE``, IID ``001``) and
  rejects IDs with several underscores. ``--double-id`` keeps the whole ID
  as both FID and IID.
* On export, ``--recode vcf`` joins FID and IID with an underscore and
  writes A2 as REF, where A2 may have been reset to the major allele.
  ``vcf-iid`` writes the IID only, and ``--keep-allele-order`` keeps the
  allele order from import, so a VCF -> PLINK -> VCF round trip preserves
  sample IDs and REF/ALT.
"""

from __future__ import annotations

import gzip
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

Runner = Callable[[list[str]], "subprocess.CompletedProcess[str]"]


def subprocess_runner(argv: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(argv, capture_output=True, text=True, check=False)


@dataclass
class ConversionResult:
    """Result of a format conversion."""

    input_path: str = ""
    output_path: str = ""
    input_format: str = ""
    output_format: str = ""
    sample_count: int = 0
    variant_count: int = 0
    success: bool = True
    error: str | None = None


def _vcf_prefix(path: str | Path) -> tuple[str, bool]:
    """Strip .vcf or .vcf.gz and say whether compressed output was asked for."""
    text = str(path)
    for suffix, compressed in ((".vcf.gz", True), (".vcf", False)):
        if text.endswith(suffix):
            return text[: -len(suffix)], compressed
    return text, False


class FormatConverter:
    """Convert between PLINK binary filesets and VCF using PLINK 1.9."""

    def __init__(self, plink_path: str = "plink", runner: Runner | None = None) -> None:
        self.plink_path = plink_path
        self.runner = runner or subprocess_runner

    def plink_to_vcf(
        self, bfile: str | Path, output_vcf: str | Path
    ) -> ConversionResult:
        """Write ``output_vcf`` (``.vcf``, or ``.vcf.gz`` for BGZF) from a fileset."""
        prefix, compressed = _vcf_prefix(output_vcf)
        written = prefix + (".vcf.gz" if compressed else ".vcf")
        result = ConversionResult(
            input_path=str(bfile),
            output_path=written,
            input_format="plink",
            output_format="vcf",
        )
        recode = ["--recode", "vcf-iid"] + (["bgz"] if compressed else [])
        argv = [
            self.plink_path,
            "--bfile",
            str(bfile),
            *recode,
            "--keep-allele-order",
            "--out",
            prefix,
        ]
        proc = self.runner(argv)
        if proc.returncode != 0:
            result.success = False
            result.error = proc.stderr or proc.stdout
            return result
        if Path(written).exists():
            result.sample_count = self.count_vcf_samples(written)
            result.variant_count = self.count_vcf_variants(written)
        return result

    def vcf_to_plink(
        self, vcf_path: str | Path, output_prefix: str | Path
    ) -> ConversionResult:
        """Write ``<output_prefix>.bed/.bim/.fam`` from a VCF."""
        result = ConversionResult(
            input_path=str(vcf_path),
            output_path=str(output_prefix),
            input_format="vcf",
            output_format="plink",
        )
        argv = [
            self.plink_path,
            "--vcf",
            str(vcf_path),
            "--double-id",
            "--keep-allele-order",
            "--make-bed",
            "--out",
            str(output_prefix),
        ]
        proc = self.runner(argv)
        if proc.returncode != 0:
            result.success = False
            result.error = proc.stderr or proc.stdout
            return result
        for suffix, attr in ((".fam", "sample_count"), (".bim", "variant_count")):
            path = Path(f"{output_prefix}{suffix}")
            if path.exists():
                with open(path) as fh:
                    setattr(result, attr, sum(1 for line in fh if line.strip()))
        return result

    @staticmethod
    def count_vcf_samples(vcf_path: str | Path) -> int:
        """Samples named on the ``#CHROM`` header line (plain or gzipped VCF)."""
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    return max(0, len(line.rstrip("\n").split("\t")) - 9)
        return 0

    @staticmethod
    def count_vcf_variants(vcf_path: str | Path) -> int:
        """Data lines in a plain or gzipped VCF."""
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt") as fh:
            return sum(1 for line in fh if line.strip() and not line.startswith("#"))
