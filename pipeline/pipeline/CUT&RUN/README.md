# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 03 | `03_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 04 | `04_call_peaks_and_define_binding.py` | Call and filter peaks, then define promoter-bound genes. |
| 05 | `05_render_figures.py` | Render and publish all registered CUT&RUN figures. |

Run this branch through the parent workflow:

```bash
python ../run_pipeline.py --step cutrun --threads 16
```

Add `--source raw` to rebuild alignment, BAM, bigWig, and peak intermediates.
The branch uses the 20200923 and 20200929 biological replicates with their
matched IgG libraries. MACS outputs are in `cutrun_work/04_peaks`, and final
promoter-gene tables are in `cutrun_work/05_promoters`.
The `figures_rendering/` directory contains required rendering modules imported
by Stage 05. They are implementation code and are not run as separate stages.
