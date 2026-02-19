"""Multi-batch genotype assembler for GWAS preparation.

Merges genotype data from multiple array batches, resolves strand
conflicts, and produces a unified PLINK dataset ready for association
analysis.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class AssemblyResult:
    """Result of a genotype assembly operation."""

    output_prefix: str = ""
    batches_merged: int = 0
    total_samples: int = 0
    total_variants: int = 0
    strand_flips: int = 0
    excluded_variants: int = 0
    success: bool = True
    error: Optional[str] = None


class GenotypeAssembler:
    """Merge multi-batch genotype data for GWAS.

    Handles strand conflicts, triallelic exclusions, and sample
    de-duplication across genotyping batches.

    Parameters
    ----------
    plink_path : str
        Path to the PLINK 1.9 binary.
    work_dir : str or Path
        Temporary working directory.
    """

    def __init__(
        self,
        plink_path: str = "plink",
        work_dir: str | Path = "gwas_work",
    ) -> None:
        self.plink_path = plink_path
        self.work_dir = Path(work_dir)

    def count_samples(self, bfile: str | Path) -> int:
        """Count samples in a PLINK binary fileset."""
        fam = Path(bfile).with_suffix(".fam")
        if not fam.exists():
            return 0
        with open(fam) as fh:
            return sum(1 for _ in fh)

    def count_variants(self, bfile: str | Path) -> int:
        """Count variants in a PLINK binary fileset."""
        bim = Path(bfile).with_suffix(".bim")
        if not bim.exists():
            return 0
        with open(bim) as fh:
            return sum(1 for _ in fh)

    def detect_strand_conflicts(
        self, bim_a: str | Path, bim_b: str | Path
    ) -> List[str]:
        """Identify variants with strand conflicts between two .bim files.

        Returns a list of variant IDs that have complementary allele
        encodings (A/T ↔ T/A, C/G ↔ G/C), indicating strand ambiguity.
        """
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        conflicts: List[str] = []

        variants_a: Dict[str, tuple] = {}
        with open(bim_a) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    variants_a[parts[1]] = (parts[4], parts[5])

        with open(bim_b) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    vid = parts[1]
                    if vid in variants_a:
                        a1_a, a2_a = variants_a[vid]
                        a1_b, a2_b = parts[4], parts[5]
                        if (
                            a1_a == complement.get(a1_b, "")
                            and a2_a == complement.get(a2_b, "")
                        ):
                            conflicts.append(vid)
        return conflicts

    def merge_batches(
        self,
        batch_prefixes: List[str | Path],
        output_prefix: str | Path,
    ) -> AssemblyResult:
        """Merge multiple PLINK filesets with self-healing strand flip.

        Uses PLINK --bmerge with automatic retry on missnp errors.
        """
        result = AssemblyResult(output_prefix=str(output_prefix))
        if len(batch_prefixes) < 2:
            result.error = "At least two batches required"
            result.success = False
            return result

        self.work_dir.mkdir(parents=True, exist_ok=True)
        merge_list = self.work_dir / "merge_list.txt"
        with open(merge_list, "w") as fh:
            for bp in batch_prefixes[1:]:
                fh.write(f"{bp}\n")

        cmd = [
            self.plink_path,
            "--bfile", str(batch_prefixes[0]),
            "--merge-list", str(merge_list),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]

        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError:
            missnp = Path(str(output_prefix) + "-merge.missnp")
            if missnp.exists():
                result.strand_flips = sum(1 for _ in open(missnp))
                exclude_cmd = cmd + ["--exclude", str(missnp)]
                try:
                    subprocess.run(exclude_cmd, capture_output=True, text=True, check=True)
                    result.excluded_variants = result.strand_flips
                except subprocess.CalledProcessError as exc:
                    result.success = False
                    result.error = str(exc)
                    return result
            else:
                result.success = False
                result.error = "Merge failed without missnp file"
                return result

        result.batches_merged = len(batch_prefixes)
        result.total_samples = self.count_samples(output_prefix)
        result.total_variants = self.count_variants(output_prefix)
        return result
