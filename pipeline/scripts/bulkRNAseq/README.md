# Bulk RNA-seq

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_prepare_inputs.py` | Prepare references, sample manifests, and local inputs. |
| 02 | `02_quantify_and_test.py` | Quantify transcripts and run differential-expression analysis. Calls Stage 03. |
| 03 | `03_render_figures.py` | Render all bulk RNA-seq figures from newly generated DEG outputs. |

## Analysis conditions

- Paired-end reads are quantified against the GENCODE M25 mouse transcriptome
  with Salmon using automatic library-type detection (`-l A`) and
  `--validateMappings`.
- Transcript-level estimates are summarized to genes with tximport using the
  generated GENCODE M25 transcript-to-gene map. Gene counts are normalized with
  edgeR trimmed mean of M-values (TMM) normalization.
- The MCM3 contrast uses the selected 20161213 libraries: non-targeting controls,
  all `sh1` libraries, and `sh2` libraries `L002` and `L003`. Genes must have
  CPM >=1 in at least two selected samples before testing.
- NONO and PSPC1 use their accepted 20180718 knockdown and non-targeting-control
  libraries. Their expression filter and TMM factors are established on the fixed
  22-library NONO/PSPC1 normalization cohort with an edgeR design-aware
  `filterByExpr` filter before target-specific contrasts are fitted.
- Differential expression uses limma linear models for knockdown versus control,
  with `voom` for MCM3 and `voomWithQualityWeights` for NONO and PSPC1. Limma
  empirical-Bayes moderation is applied to the fitted contrasts.
- P-values are Benjamini-Hochberg adjusted. Differentially expressed genes require
  **BH-adjusted p-value <= 0.05** and **absolute log2 fold change >= 0.28**.

## Outputs

Stage 03 creates per-target DEG tables and the canonical `Volcano_*`, `Heatmap_*`,
`VennDiagram_UP*`, and `VennDiagram_DOWN*` PNG files under `deg_work`. The
`*_noTexts` Venn variants intentionally omit all labels and values inside the
diagram.

Run Stage 03 with `--publication-figures` to render text-free versions of every
panel in `deg_work/visuals/publication_figures`.
