"""Format conversion module.

Converts between PLINK binary (.bed/.bim/.fam) and VCF formats
using PLINK and bcftools.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


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
    error: Optional[str] = None


class FormatConverter:
    """Convert between PLINK and VCF formats.

    Parameters
    ----------
    plink_path : str
        Path to PLINK 1.9 binary.
    bcftools_path : str
        Path to bcftools binary.
    """

    def __init__(
        self,
        plink_path: str = "plink",
        bcftools_path: str = "bcftools",
    ) -> None:
        self.plink_path = plink_path
        self.bcftools_path = bcftools_path

    def plink_to_vcf(
        self,
        bfile: str | Path,
        output_vcf: str | Path,
    ) -> ConversionResult:
        """Convert PLINK binary fileset to VCF."""
        result = ConversionResult(
            input_path=str(bfile),
            output_path=str(output_vcf),
            input_format="plink",
            output_format="vcf",
        )
        cmd = [
            self.plink_path,
            "--bfile", str(bfile),
            "--recode", "vcf",
            "--out", str(Path(output_vcf).with_suffix("")),
            "--allow-no-sex",
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            result.success = True
        except subprocess.CalledProcessError as exc:
            result.success = False
            result.error = exc.stderr
        return result

    def vcf_to_plink(
        self,
        vcf_path: str | Path,
        output_prefix: str | Path,
    ) -> ConversionResult:
        """Convert VCF to PLINK binary fileset."""
        result = ConversionResult(
            input_path=str(vcf_path),
            output_path=str(output_prefix),
            input_format="vcf",
            output_format="plink",
        )
        cmd = [
            self.plink_path,
            "--vcf", str(vcf_path),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            result.success = True
        except subprocess.CalledProcessError as exc:
            result.success = False
            result.error = exc.stderr
        return result

    @staticmethod
    def count_vcf_samples(vcf_path: str | Path) -> int:
        """Count samples in a VCF file from the header line."""
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    parts = line.strip().split("\t")
                    return max(0, len(parts) - 9)
        return 0

    @staticmethod
    def count_vcf_variants(vcf_path: str | Path) -> int:
        """Count non-header lines in a VCF file."""
        count = 0
        with open(vcf_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    count += 1
        return count
