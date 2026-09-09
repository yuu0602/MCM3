# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 03 | `03_cutrun.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 04 | `04_promoter_binding.py` | Call and filter peaks, then define promoter-bound genes. |
| 05 | `05_cutrun_figures.py` | Render and publish all registered CUT&RUN figures. |

Run this branch through the parent workflow:

```bash
python ../run_all.py --step cutrun --threads 16
```

Add `--source raw` to rebuild alignment, BAM, bigWig, and peak intermediates.
MACS outputs from 2019 and 2020 share `cutrun_work/04_peaks`, separated only by
year. Final promoter-gene tables are in `cutrun_work/05_promoters`.
The `figures_rendering/` directory contains required rendering modules imported
by Stage 05. They are implementation code and are not run as separate stages.
