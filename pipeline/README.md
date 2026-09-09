# Pipeline

`run_all.py` executes the numbered stages in dependency order.

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare.py` | References, manifests, and local inputs |
| 02 | `bulkRNAseq/02_rnaseq.py` | DEG tables and RNA-seq figures |
| 03 | `CUT&RUN/03_cutrun.py` | Alignments, filtered BAMs, and normalized tracks |
| 04 | `CUT&RUN/04_promoter_binding.py` | Peaks and promoter-bound genes |
| 05 | `CUT&RUN/05_cutrun_figures.py` | CUT&RUN figures |
| 06 | `regulatory_target/06_regulatory_targets.py` | Regulatory targets and figures |
| 07 | `regulatory_target/07_verify.py` | Output and consistency checks |

Run all stages:

```bash
python run_all.py --source accepted --step all --threads 16
```

Use `--source raw` for a complete FASTQ rebuild. Figure helpers are kept in each
assay's `figures_rendering/` directory and are called by the numbered parent
script.
