"""Merge genotyping batches into one PLINK fileset, with strand checks.

Each later batch is compared with the first (the reference) variant by
variant, using the alleles in the .bim files:

* **flipped**: a non-ambiguous SNP whose alleles match only after taking the
  complement (A<->T, C<->G). Fixed with ``plink --flip`` on that batch.
* **ambiguous**: an A/T or C/G SNP. Its strand cannot be told from the
  alleles, so it is reported and, if requested, excluded.
* **incompatible**: alleles that do not match even after flipping (for
  example A/G against A/C). Always excluded.

Only then are the batches merged, once, with ``plink --merge-list``. This
avoids PLINK's trial-and-error loop of flipping every variant in the
``.missnp`` file, which with several batches can flip a variant in a batch
that was already correct.
"""

from __future__ import annotations

import subprocess
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

Runner = Callable[[list[str]], int]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def subprocess_runner(argv: list[str]) -> int:
    """Run a command quietly and return its exit code."""
    return subprocess.run(argv, capture_output=True, text=True, check=False).returncode


def _read_bim(path: str | Path) -> dict[str, tuple[str, str]]:
    alleles: dict[str, tuple[str, str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 6:
                alleles[parts[1]] = (parts[4].upper(), parts[5].upper())
    return alleles


def _count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


@dataclass
class StrandCheck:
    """Strand comparison of one batch against the reference batch."""

    flipped: list[str] = field(default_factory=list)
    ambiguous: list[str] = field(default_factory=list)
    incompatible: list[str] = field(default_factory=list)


def compare_alleles(ref: tuple[str, str], other: tuple[str, str]) -> str:
    """Classify one shared variant: 'same', 'flipped', 'ambiguous' or 'incompatible'.

    Allele order is ignored (PLINK may store A1/A2 either way round), and a
    missing allele (``0``, for a variant monomorphic in that batch) matches
    whatever the other batch has.
    """
    ref_set = {a for a in ref if a != "0"}
    other_set = {a for a in other if a != "0"}
    union = ref_set | other_set
    if union in ({"A", "T"}, {"C", "G"}):
        return "ambiguous"
    if len(union) <= 2:
        return "same"
    flipped = {COMPLEMENT.get(a, "?") for a in other_set}
    if "?" not in flipped and len(ref_set | flipped) <= 2:
        return "flipped"
    return "incompatible"


@dataclass
class AssemblyResult:
    """Outcome of a merge."""

    output_prefix: str = ""
    batches_merged: int = 0
    total_samples: int = 0
    total_variants: int = 0
    strand_flips: int = 0
    ambiguous_variants: int = 0
    excluded_variants: int = 0
    success: bool = True
    error: str | None = None
    commands: list[list[str]] = field(default_factory=list)


class GenotypeAssembler:
    """Merge PLINK binary filesets from several genotyping batches.

    Parameters
    ----------
    plink_path : str
        PLINK 1.9 executable.
    work_dir : path
        Directory for flipped and filtered intermediate filesets.
    runner : callable, optional
        Runs one command and returns its exit code (for testing).
    """

    def __init__(
        self,
        plink_path: str = "plink",
        work_dir: str | Path = "gwas_work",
        runner: Runner | None = None,
    ) -> None:
        self.plink_path = plink_path
        self.work_dir = Path(work_dir)
        self.runner = runner or subprocess_runner

    @staticmethod
    def count_samples(bfile: str | Path) -> int:
        """Samples in ``<bfile>.fam`` (0 if the file is missing)."""
        return _count_lines(Path(f"{bfile}.fam"))

    @staticmethod
    def count_variants(bfile: str | Path) -> int:
        """Variants in ``<bfile>.bim`` (0 if the file is missing)."""
        return _count_lines(Path(f"{bfile}.bim"))

    @staticmethod
    def check_strand(bim_ref: str | Path, bim_other: str | Path) -> StrandCheck:
        """Compare shared variants between two .bim files."""
        ref, other = _read_bim(bim_ref), _read_bim(bim_other)
        check = StrandCheck()
        for vid in sorted(ref.keys() & other.keys()):
            kind = compare_alleles(ref[vid], other[vid])
            if kind != "same":
                getattr(check, kind).append(vid)
        return check

    def detect_strand_conflicts(
        self, bim_a: str | Path, bim_b: str | Path
    ) -> list[str]:
        """Shared variants whose strand cannot be trusted as is.

        Flipped SNPs plus A/T and C/G SNPs, sorted.
        """
        check = self.check_strand(bim_a, bim_b)
        return sorted(check.flipped + check.ambiguous)

    def _run(self, result: AssemblyResult, argv: list[str]) -> bool:
        result.commands.append(argv)
        return self.runner(argv) == 0

    def merge_batches(
        self,
        batch_prefixes: Sequence[str | Path],
        output_prefix: str | Path,
        exclude_ambiguous: bool = False,
    ) -> AssemblyResult:
        """Flip, filter and merge batches into ``output_prefix``.

        Parameters
        ----------
        batch_prefixes : sequence of paths
            PLINK fileset prefixes; the first is the strand reference.
        output_prefix : path
            Prefix for the merged fileset.
        exclude_ambiguous : bool
            Also drop A/T and C/G SNPs. Advisable when batches come from
            different arrays or pipelines; usually unnecessary when they
            share one array and one strand convention.
        """
        result = AssemblyResult(output_prefix=str(output_prefix))
        if len(batch_prefixes) < 2:
            result.success = False
            result.error = "At least two batches required"
            return result

        self.work_dir.mkdir(parents=True, exist_ok=True)
        ref = str(batch_prefixes[0])
        exclude: set[str] = set()
        ambiguous: set[str] = set()
        prepared = [ref]
        for i, prefix in enumerate(batch_prefixes[1:], start=1):
            prefix = str(prefix)
            check = self.check_strand(f"{ref}.bim", f"{prefix}.bim")
            ambiguous.update(check.ambiguous)
            exclude.update(check.incompatible)
            if check.flipped:
                flip_list = self.work_dir / f"batch{i}.flip"
                flip_list.write_text("\n".join(check.flipped) + "\n")
                flipped_prefix = str(self.work_dir / f"batch{i}_flipped")
                argv = [
                    self.plink_path,
                    "--bfile",
                    prefix,
                    "--flip",
                    str(flip_list),
                    "--make-bed",
                    "--out",
                    flipped_prefix,
                ]
                if not self._run(result, argv):
                    result.success = False
                    result.error = f"plink --flip failed for {prefix}"
                    return result
                result.strand_flips += len(check.flipped)
                prefix = flipped_prefix
            prepared.append(prefix)

        result.ambiguous_variants = len(ambiguous)
        if exclude_ambiguous:
            exclude |= ambiguous
        if exclude:
            exclude_list = self.work_dir / "exclude.txt"
            exclude_list.write_text("\n".join(sorted(exclude)) + "\n")
            filtered = []
            for i, prefix in enumerate(prepared):
                out = str(self.work_dir / f"batch{i}_filtered")
                argv = [
                    self.plink_path,
                    "--bfile",
                    prefix,
                    "--exclude",
                    str(exclude_list),
                    "--make-bed",
                    "--out",
                    out,
                ]
                if not self._run(result, argv):
                    result.success = False
                    result.error = f"plink --exclude failed for {prefix}"
                    return result
                filtered.append(out)
            prepared = filtered
            result.excluded_variants = len(exclude)

        merge_list = self.work_dir / "merge_list.txt"
        merge_list.write_text("\n".join(prepared[1:]) + "\n")
        argv = [
            self.plink_path,
            "--bfile",
            prepared[0],
            "--merge-list",
            str(merge_list),
            "--make-bed",
            "--out",
            str(output_prefix),
        ]
        if not self._run(result, argv):
            missnp = Path(f"{output_prefix}-merge.missnp")
            result.success = False
            result.error = (
                f"merge failed; {_count_lines(missnp)} variants listed in {missnp}"
                if missnp.exists()
                else "merge failed without a .missnp file"
            )
            return result

        result.batches_merged = len(batch_prefixes)
        result.total_samples = self.count_samples(output_prefix)
        result.total_variants = self.count_variants(output_prefix)
        return result
