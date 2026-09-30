"""GWAS quality control from PLINK 1.9 report files.

Each check reads one PLINK output and returns the IDs that fail:

========================  ==========================  =====================
Check                     PLINK command               File
========================  ==========================  =====================
Variant call rate         ``--missing``               ``.lmiss``
Sample call rate          ``--missing``               ``.imiss``
Minor allele frequency    ``--freq``                  ``.frq``
Hardy-Weinberg            ``--hardy``                 ``.hwe``
Heterozygosity outliers   ``--het``                   ``.het``
Sex discordance           ``--check-sex``             ``.sexcheck``
Relatedness               ``--genome``                ``.genome``
========================  ==========================  =====================

Columns are found by header name, so PLINK's padded, whitespace-aligned
output is read correctly.
"""

from __future__ import annotations

import statistics
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path


def _rows(path: str | Path) -> Iterator[dict[str, str]]:
    """Yield each data row of a PLINK report as ``{column: value}``."""
    with open(path) as fh:
        header = fh.readline().split()
        for line in fh:
            parts = line.split()
            if len(parts) == len(header):
                yield dict(zip(header, parts, strict=True))


def _float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None  # "NA", "nan" or similar


@dataclass
class QCReport:
    """Counts before and after QC, and removals per filter."""

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
    """Threshold PLINK QC reports.

    Parameters
    ----------
    call_rate_variant, call_rate_sample : float
        Minimum call rates (default 0.98).
    maf_threshold : float
        Minimum minor allele frequency (default 0.01).
    hwe_p : float
        Hardy-Weinberg exact test p-value below which a variant fails
        (default 1e-6).
    het_sd : float
        Standard deviations from the mean inbreeding coefficient F beyond
        which a sample fails (default 3).
    pi_hat : float
        Relatedness threshold (default 0.2, roughly second-degree relatives).
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

    def check_variant_call_rates(self, lmiss_path: str | Path) -> set[str]:
        """Variants whose call rate (1 - F_MISS) is below the threshold."""
        return {
            r["SNP"]
            for r in _rows(lmiss_path)
            if (f := _float(r["F_MISS"])) is not None
            and 1.0 - f < self.call_rate_variant
        }

    def check_sample_call_rates(self, imiss_path: str | Path) -> set[str]:
        """Samples (IID) whose call rate is below the threshold."""
        return {
            r["IID"]
            for r in _rows(imiss_path)
            if (f := _float(r["F_MISS"])) is not None
            and 1.0 - f < self.call_rate_sample
        }

    def check_maf(self, frq_path: str | Path) -> set[str]:
        """Variants with MAF below the threshold (monomorphic 'NA' included)."""
        failed: set[str] = set()
        for r in _rows(frq_path):
            maf = _float(r["MAF"])
            if maf is None or maf < self.maf_threshold:
                failed.add(r["SNP"])
        return failed

    def check_hwe(self, hwe_path: str | Path, test: str = "ALL") -> set[str]:
        """Variants failing Hardy-Weinberg equilibrium.

        ``test`` picks the row type PLINK reports: ``ALL`` for all samples,
        or ``UNAFF`` to test controls only in a case-control study, which
        avoids removing true associations that distort HWE in cases.
        """
        return {
            r["SNP"]
            for r in _rows(hwe_path)
            if r["TEST"] == test
            and (p := _float(r["P"])) is not None
            and p < self.hwe_p
        }

    def check_heterozygosity(self, het_path: str | Path) -> set[str]:
        """Samples whose inbreeding coefficient F is more than het_sd SDs from the mean.

        Very low F (excess heterozygosity) often means contamination; very
        high F can mean poor DNA quality or consanguinity.
        """
        f_values = {
            r["IID"]: f for r in _rows(het_path) if (f := _float(r["F"])) is not None
        }
        if len(f_values) < 2:
            return set()
        mean = statistics.fmean(f_values.values())
        sd = statistics.stdev(f_values.values())
        return {iid for iid, f in f_values.items() if abs(f - mean) > self.het_sd * sd}

    @staticmethod
    def check_sex(sexcheck_path: str | Path) -> set[str]:
        """Samples PLINK marks as ``PROBLEM`` (reported and genetic sex differ)."""
        return {r["IID"] for r in _rows(sexcheck_path) if r["STATUS"] == "PROBLEM"}

    def check_relatedness(self, genome_path: str | Path) -> set[str]:
        """Samples to remove so that no pair has PI_HAT above the threshold.

        Greedy: repeatedly remove the sample involved in the most remaining
        related pairs (ties broken by ID), which usually removes fewer people
        than dropping one member of every pair at random.
        """
        pairs: set[tuple[str, str]] = set()
        for r in _rows(genome_path):
            pi = _float(r["PI_HAT"])
            if pi is not None and pi > self.pi_hat:
                pairs.add((r["IID1"], r["IID2"]))
        removed: set[str] = set()
        while pairs:
            degree: dict[str, int] = {}
            for a, b in pairs:
                degree[a] = degree.get(a, 0) + 1
                degree[b] = degree.get(b, 0) + 1
            worst = min(degree, key=lambda s: (-degree[s], s))
            removed.add(worst)
            pairs = {p for p in pairs if worst not in p}
        return removed

    def generate_report(
        self,
        initial_samples: int,
        initial_variants: int,
        failed_variants: dict[str, set[str]],
        failed_samples: dict[str, set[str]],
    ) -> QCReport:
        """Summarise removals. Keys: call_rate, maf, hwe for variants;
        call_rate, heterozygosity, sex, relatedness for samples. Final counts
        subtract the union, so an ID failing two filters is removed once."""
        fv = {k: failed_variants.get(k, set()) for k in ("call_rate", "maf", "hwe")}
        fs = {
            k: failed_samples.get(k, set())
            for k in ("call_rate", "heterozygosity", "sex", "relatedness")
        }
        return QCReport(
            initial_samples=initial_samples,
            initial_variants=initial_variants,
            removed_low_call_rate_variants=len(fv["call_rate"]),
            removed_low_maf_variants=len(fv["maf"]),
            removed_hwe_variants=len(fv["hwe"]),
            removed_low_call_rate_samples=len(fs["call_rate"]),
            removed_het_outlier_samples=len(fs["heterozygosity"]),
            removed_sex_discordance_samples=len(fs["sex"]),
            removed_related_samples=len(fs["relatedness"]),
            final_variants=initial_variants - len(set().union(*fv.values())),
            final_samples=initial_samples - len(set().union(*fs.values())),
        )
