# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. |
| 03 | `03_render_figures.py` | Render figures; use `--igv-tracks` to also generate RNA-seq IGV tracks. |
