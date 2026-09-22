# MCM3 Project Pipeline

This repository contains the publication pipeline for the MCM3-PSPC1-NONO project.

## Workflow

1. `bulkRNAseq/01_prepare_inputs.py`: bulk-RNAseq input and reference preparation.
2. `bulkRNAseq/02_quantify_and_test.py`: bulk-RNAseq quantification and differential expression.
3. `bulkRNAseq/03_render_figures.py`: bulk-RNAseq figure rendering.
4. `CUT&RUN/04_align_and_normalize.py`: alignment, filtering, yeast normalization, and coverage tracks.
5. `CUT&RUN/05_call_peaks_and_define_binding.py`: matched-IgG peak calling, filtering, and promoter binding.
6. `CUT&RUN/06_render_figures.py`: CUT&RUN figure rendering.
7. `regulatory_target/07_define_regulatory_targets.py`: regulatory target analysis and figure rendering.
