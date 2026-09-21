# MCM3 Project Analysis Pipeline
## Analysis Conditions

- RNA-seq reads are quantified with Salmon and tested with limma-voom. DEGs
  require **BH-adjusted p-value <= 0.05** and **absolute log2 fold change >= 0.28**.
- CUT&RUN retains properly paired, primary, non-duplicate mouse alignments with
  **MAPQ >= 30**. Coverage tracks are normalized with retained yeast fragments.
- MACS3 candidate regions are called at **p <= 1e-3**. Pooled matched-IgG peaks
  retained downstream require **p <= 1e-4, q <= 0.01, and fold enrichment >= 3**.
- Promoters are **GENCODE M25 TSS +/- 1,000 bp**. A peak must overlap a promoter by
  **at least 250 bp**. Final binding requires a **promoter factor/IgG ratio >= 2** in
  both biological replicates.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.
