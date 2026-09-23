# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. |
| 03 | `03_render_figures.py` | Render figures; use `--igv-tracks` to also generate RNA-seq IGV tracks. |

`03_render_figures.py --igv-tracks` uses the selected libraries in
`metadata/Samples.tsv` to generate CPM-normalized, 10-bp individual bigWigs and
equal-weight `NT_mean.bw`, `shMCM3_mean.bw`, `shNONO_mean.bw`, and
`shPSPC1_mean.bw` tracks in `deg_work/igv_tracks/`.
