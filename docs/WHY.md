# Why it's built this way

## The problem

GWAS preparation is a sequence of routine steps (merge batches, filter variants and samples, convert formats), each with a known way of going quietly wrong. The aim is to make each step explicit, testable and honest about what it did.

## Design choices

**Why compare each batch with a reference batch instead of retrying PLINK merges?** Because the `.missnp` file only says which variants had more than two alleles across the merge, not which batch is wrong. Flipping them everywhere works for two datasets, but with several batches it can flip a SNP in a batch that was already right. Comparing alleles batch by batch says exactly where to flip.

**Why treat A/T and C/G SNPs separately?** Because complementing A/T gives T/A: the alleles cannot show whether the strand differs. Within one array and pipeline they are usually safe; across arrays they are not. The code reports them and leaves the exclusion decision to the analyst.

**Why read PLINK reports by column name?** Because PLINK pads columns with spaces and different commands put fields in different places. Looking up `F_MISS`, `MAF`, `P`, `F`, `STATUS` or `PI_HAT` by name survives both.

**Why greedy removal for relatedness?** Because one person related to three others should be removed once, not three others removed. Taking out the most-connected sample first usually keeps more of the cohort.

**Why offer HWE on controls only?** Because a true disease association can distort genotype frequencies in cases. Testing HWE across everyone can remove the very signal the study is looking for.

**Why `--double-id`, `vcf-iid` and `--keep-allele-order` in conversions?** Because PLINK 1.9's defaults split VCF sample IDs at an underscore on import, join FID and IID with one on export, and may write the major allele as REF. Each default is reasonable on its own; together they make a round trip lose information.

**Why an injectable runner for PLINK?** Because the logic worth testing is the plan (which batch is flipped, what is excluded, what happens when the merge fails), and that can be tested without PLINK or real data.

## Questions worth asking

**"Why not resolve ambiguous SNPs using allele frequencies?"**
That is the usual next step: if an A/T SNP has frequency 0.2 in one batch and 0.8 in the other, it is probably flipped. It is not implemented yet because it needs a frequency threshold and breaks down near 0.5, where both orientations look alike. It is the first item on the roadmap.

**"Is PI_HAT above 0.2 the right relatedness cut-off?"**
It is a common default that removes second-degree relatives and closer. The right threshold depends on the analysis: mixed models can keep relatives, while standard logistic regression cannot. The threshold is a parameter for that reason.

**"You read PLINK's output. What if someone ran PLINK with different settings?"**
The checks only threshold what PLINK reported, so they inherit its settings (for example LD pruning before `--het` and `--genome`, which should be done). The functions document which PLINK command each file comes from; recording the exact commands alongside the reports would make that auditable.

## What's next

- Resolve ambiguous SNPs by allele frequency against the reference batch.
- Add principal-component ancestry checks.
- Record the PLINK commands and versions alongside each QC run.
