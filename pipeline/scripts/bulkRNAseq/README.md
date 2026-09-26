# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. |
| 03 | `03_render_figures.py` | Render figures; `--igv-tracks` additionally uses HISAT2 |

Because bigwig file generation for bulkRNAseq data is memory-intensive, regular laptops (e.g. 24GB memory) cannot handle the computational process. We suggest you use HPC for this; otherwize you would skip it.
