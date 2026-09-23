# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. |
| 03 | `03_render_figures.py` | Render all bulk RNA-seq figures from newly generated DEG outputs. |
