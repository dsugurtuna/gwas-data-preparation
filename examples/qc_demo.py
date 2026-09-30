"""Run every QC check and a strand comparison on the synthetic example files.

python examples/qc_demo.py
"""

from pathlib import Path

from gwas_prep import GenotypeAssembler, QualityController

HERE = Path(__file__).parent

qc = QualityController()
failed_variants = {
    "call_rate": qc.check_variant_call_rates(HERE / "cohort.lmiss"),
    "maf": qc.check_maf(HERE / "cohort.frq"),
    "hwe": qc.check_hwe(HERE / "cohort.hwe"),
}
failed_samples = {
    "call_rate": qc.check_sample_call_rates(HERE / "cohort.imiss"),
    "heterozygosity": qc.check_heterozygosity(HERE / "cohort.het"),
    "sex": qc.check_sex(HERE / "cohort.sexcheck"),
    "relatedness": qc.check_relatedness(HERE / "cohort.genome"),
}
for kind, failed in (("variant", failed_variants), ("sample", failed_samples)):
    for name, ids in failed.items():
        print(f"{kind} {name:15s} {sorted(ids)}")

report = qc.generate_report(40, 60, failed_variants, failed_samples)
print(
    f"variants kept {report.final_variants}/{report.initial_variants}, "
    f"samples kept {report.final_samples}/{report.initial_samples}"
)

check = GenotypeAssembler.check_strand(HERE / "batch1.bim", HERE / "batch2.bim")
print(
    f"strand: flipped {check.flipped}, ambiguous {check.ambiguous}, "
    f"incompatible {check.incompatible}"
)
