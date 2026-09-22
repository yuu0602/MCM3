# MCM3 Project Analysis Pipeline
## Analysis Conditions

- bulk-RNAseq reads are quantified with Salmon and tested with limma-voom. DEGs
  require **BH-adjusted p-value <= 0.05** and **absolute log2 fold change >= 0.28**.
- CUT&RUN mouse and yeast reads are aligned separately with Bowtie2. Properly
  paired primary mouse alignments with **MAPQ >= 30** are retained.
- Mouse paired-fragment coverage is yeast-normalized to **10,000 retained yeast
  read pairs** per library. Factor and IgG libraries use the same alignment,
  filtering, and normalization procedure.
- For each factor, two biological-replicate BAMPE libraries are pooled and compared
  with pooled matched-IgG libraries using MACS3. Downstream peaks require a direct
  **q-value <= 0.05** and **fold enrichment >= 3**; no separate p-value threshold is
  applied.
- Promoters are **GENCODE M25 TSS +/- 1,000 bp**. A retained pooled peak must overlap
  a promoter by **at least 250 bp**. A promoter-bound gene additionally requires a
  yeast-normalized factor/IgG coverage ratio **>= 2** in **both biological replicates**.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.
