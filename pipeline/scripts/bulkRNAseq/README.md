# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. |
| 03 | `03_render_figures.py` | Render figures; `--igv-tracks` additionally uses HISAT2 with GENCODE M25 splice/exon guidance to generate MAPQ >=30, properly paired CPM bigWig tracks. The one-time HISAT2 index is built in the local temporary filesystem and copied to `deg_work/igv_tracks/reference`; set `MCM3_HISAT2_INDEX_TMPDIR` to override that workspace. |
