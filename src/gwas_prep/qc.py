"""Quality control module for GWAS datasets.

Applies standard GWAS QC filters: call rate, MAF, HWE, heterozygosity,
sex discordance, and relatedness.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Set


@dataclass
class QCReport:
    """Results of quality control filtering."""

    initial_samples: int = 0
    initial_variants: int = 0
    final_samples: int = 0
    final_variants: int = 0
    removed_low_call_rate_variants: int = 0
    removed_low_maf_variants: int = 0
    removed_hwe_variants: int = 0
    removed_low_call_rate_samples: int = 0
    removed_het_outlier_samples: int = 0
    removed_sex_discordance_samples: int = 0
    removed_related_samples: int = 0

    @property
    def variant_pass_rate(self) -> float:
        if self.initial_variants == 0:
            return 0.0
        return self.final_variants / self.initial_variants

    @property
    def sample_pass_rate(self) -> float:
        if self.initial_samples == 0:
            return 0.0
        return self.final_samples / self.initial_samples


class QualityController:
    """GWAS QC filter pipeline.

    Applies sequential variant and sample-level filters with
    configurable thresholds.

    Parameters
    ----------
    call_rate_variant : float
        Minimum variant call rate (default 0.98).
    call_rate_sample : float
        Minimum sample call rate (default 0.98).
    maf_threshold : float
        Minimum minor allele frequency (default 0.01).
    hwe_p : float
        Hardy-Weinberg equilibrium p-value threshold (default 1e-6).
    het_sd : float
        Heterozygosity outlier threshold in standard deviations (default 3.0).
    pi_hat : float
        Relatedness threshold (default 0.2).
    """

    def __init__(
        self,
        call_rate_variant: float = 0.98,
        call_rate_sample: float = 0.98,
        maf_threshold: float = 0.01,
        hwe_p: float = 1e-6,
        het_sd: float = 3.0,
        pi_hat: float = 0.2,
    ) -> None:
        self.call_rate_variant = call_rate_variant
        self.call_rate_sample = call_rate_sample
        self.maf_threshold = maf_threshold
        self.hwe_p = hwe_p
        self.het_sd = het_sd
        self.pi_hat = pi_hat

    def check_variant_call_rates(
        self, lmiss_path: str | Path
    ) -> Set[str]:
        """Parse PLINK .lmiss file and return variants below threshold.

        The .lmiss file columns: CHR, SNP, N_MISS, N_GENO, F_MISS
        """
        failed: Set[str] = set()
        with open(lmiss_path) as fh:
            header = True
            for line in fh:
                if header:
                    header = False
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    f_miss = float(parts[4])
                    if (1.0 - f_miss) < self.call_rate_variant:
                        failed.add(parts[1])
        return failed

    def check_sample_call_rates(
        self, imiss_path: str | Path
    ) -> Set[str]:
        """Parse PLINK .imiss file and return samples below threshold.

        The .imiss file columns: FID, IID, MISS_PHENO, N_MISS, N_GENO, F_MISS
        """
        failed: Set[str] = set()
        with open(imiss_path) as fh:
            header = True
            for line in fh:
                if header:
                    header = False
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    f_miss = float(parts[5])
                    if (1.0 - f_miss) < self.call_rate_sample:
                        failed.add(parts[1])
        return failed

    def check_maf(self, frq_path: str | Path) -> Set[str]:
        """Parse PLINK .frq file and return variants below MAF threshold.

        Columns: CHR, SNP, A1, A2, MAF, NCHROBS
        """
        failed: Set[str] = set()
        with open(frq_path) as fh:
            header = True
            for line in fh:
                if header:
                    header = False
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    maf = float(parts[4])
                    if maf < self.maf_threshold:
                        failed.add(parts[1])
        return failed

    def generate_report(
        self,
        initial_samples: int,
        initial_variants: int,
        failed_variants: Dict[str, Set[str]],
        failed_samples: Dict[str, Set[str]],
    ) -> QCReport:
        """Summarise QC filtering into a report."""
        report = QCReport(
            initial_samples=initial_samples,
            initial_variants=initial_variants,
        )

        all_failed_variants: Set[str] = set()
        report.removed_low_call_rate_variants = len(failed_variants.get("call_rate", set()))
        all_failed_variants |= failed_variants.get("call_rate", set())
        report.removed_low_maf_variants = len(failed_variants.get("maf", set()))
        all_failed_variants |= failed_variants.get("maf", set())
        report.removed_hwe_variants = len(failed_variants.get("hwe", set()))
        all_failed_variants |= failed_variants.get("hwe", set())

        all_failed_samples: Set[str] = set()
        report.removed_low_call_rate_samples = len(failed_samples.get("call_rate", set()))
        all_failed_samples |= failed_samples.get("call_rate", set())
        report.removed_het_outlier_samples = len(failed_samples.get("heterozygosity", set()))
        all_failed_samples |= failed_samples.get("heterozygosity", set())
        report.removed_sex_discordance_samples = len(failed_samples.get("sex", set()))
        all_failed_samples |= failed_samples.get("sex", set())
        report.removed_related_samples = len(failed_samples.get("relatedness", set()))
        all_failed_samples |= failed_samples.get("relatedness", set())

        report.final_variants = initial_variants - len(all_failed_variants)
        report.final_samples = initial_samples - len(all_failed_samples)
        return report
