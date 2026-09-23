# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 04 | `04_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 05 | `05_call_peaks_and_define_binding.py` | Call and filter peaks, then define promoter-bound genes. |
| 06 | `06_render_figures.py` | Render all CUT&RUN figures. |

Run Stage 06 with `--publication-figures` to render text-free versions of every
panel in `cutrun_work/visuals/publication_figures`. They use the `_noTexts` suffix.
