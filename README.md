# MCM3 Project Pipelines

This repository contains the two condition-specific publication pipelines for the MCM3-PSPC1-NONO project. The older root-level pipeline has been removed.

| Directory | Promoter peak overlap | Factor/IgG support |
| --- | --- | --- |
| `strict_2p0fold_250bp/` | >=250 bp | >=2.0 in both biological replicates |
| `relaxed_1p5fold_200bp/` | >=200 bp | >=1.5 in both biological replicates |

Each condition directory contains its own `pipeline/`, `README.md`, `environment.yml`, and `tests/test_runtime_hygiene.py`. The two pipelines use the same seven-stage layout:

1. RNA-seq input/reference preparation.
2. RNA-seq quantification, differential expression, and figures.
3. CUT&RUN alignment, filtering, yeast normalization, and tracks.
4. Matched-IgG peak calling/filtering and promoter binding.
5. CUT&RUN figures.
6. Regulatory-target intersection and figures.
7. Output validation.

Raw FASTQs, alignments, bigWigs, peak outputs, and publication figures are intentionally not tracked in GitHub. To rebuild from FASTQ files, place `bulkRNAseq_data/` and `CUT&RUN_data/` beside the two condition folders and follow the relevant condition README.

The reference transcript FASTA is downloaded by the pipeline during a raw rebuild and is also not tracked.
