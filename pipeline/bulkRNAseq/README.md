# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_rnaseq.py` | Quantify transcripts, test differential expression, and render RNA-seq figures. |

Run through the parent workflow:

```bash
python ../run_all.py --step rna --threads 16
```

Add `--source raw` to regenerate Salmon quantifications from FASTQs.
Quantifications are written to `deg_work/salmon`; native Salmon output names
are retained, while project-level tables use the corresponding figure names.

The `figures_rendering/` directory contains the required Venn renderer imported
by Stage 02. It is implementation code and is not run as a separate stage.
