# Pipeline

`run_pipeline.py` executes the numbered stages in dependency order.

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare_inputs.py` | References, manifests, and local inputs |
| 02 | `bulkRNAseq/02_quantify_and_test.py` | DEG tables and RNA-seq figures |
| 03 | `CUT&RUN/03_align_and_normalize.py` | Alignments, filtered BAMs, and normalized tracks |
| 04 | `CUT&RUN/04_call_peaks_and_define_binding.py` | Peaks and promoter-bound genes |
| 05 | `CUT&RUN/05_render_figures.py` | CUT&RUN figures |
| 06 | `regulatory_target/06_define_regulatory_targets.py` | Regulatory targets and figures |
