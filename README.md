# MCM3 Project Pipeline

This repository contains the publication pipeline for the MCM3-PSPC1-NONO project.

The retained `pipeline/` uses a promoter-peak overlap threshold of >=250 bp and requires factor/IgG support >=2.0 in both biological replicates. It contains its own `README.md`, `environment.yml`, and `tests/test_runtime_hygiene.py`, with the following seven-stage layout:

1. RNA-seq input/reference preparation.
2. RNA-seq quantification, differential expression, and figures.
3. CUT&RUN alignment, filtering, yeast normalization, and tracks.
4. Matched-IgG peak calling/filtering and promoter binding.
5. CUT&RUN figures.
6. Regulatory-target intersection and figures.
7. Output validation.

Raw FASTQs, alignments, bigWigs, peak outputs, and publication figures are intentionally not tracked in GitHub. To rebuild from FASTQ files, place `bulkRNAseq_data/` and `CUT&RUN_data/` beside `pipeline/` and follow its README.

The reference transcript FASTA is downloaded by the pipeline during a raw rebuild and is also not tracked.
