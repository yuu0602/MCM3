# MCM3 Project Analysis Pipeline

This directory contains the RNA-seq, CUT&RUN, and regulatory-target analyses
used for the final figures. All final figures are generated inside this
directory; no pipeline stage copies a figure from another reanalysis folder.

## Run The Workflow

Before creating the software environment, replace each `=VERSION` placeholder in
`environment.yml` with a manuscript-verified version. No manuscript was available
to establish these versions; the placeholder file is not yet installable.
Existing pins are retained, including the minor-series pins for Python and R.

Create the software environment:

```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml
conda activate mcm3-pipeline
```

On Apple Silicon, this uses the compatible Intel Conda package set through
Rosetta; the pinned Salmon and Bioconductor dependencies cannot be solved
together from the native ARM channels.

Rebuild tables and figures from the local accepted alignment, peak, and Salmon
intermediates:

```bash
python pipeline/run_pipeline.py --source accepted --step all --threads 16
```

Rebuild from raw FASTQs, including RNA-seq quantification and CUT&RUN
alignment, filtering, track generation, and peak calling:

```bash
python pipeline/run_pipeline.py --source raw --step all --threads 16
```

Raw FASTQs are stored beside this directory in the `bulkRNAseq_data/` and
`CUT&RUN_data/` directories. Discovery is recursive, so date and sample
subdirectories are supported; filenames are underscore-delimited. CUT&RUN
uses only the 20200923 and 20200929 samples.

Use `--dry-run` to print the stages without running them.

## Analysis Summary

- RNA-seq reads are quantified with Salmon and tested with limma-voom. DEGs
  require BH-adjusted p-value <= 0.05 and absolute log2 fold change >= 0.28.
- CUT&RUN retains properly paired, primary, non-duplicate mouse alignments with
  MAPQ >= 30. Coverage tracks are normalized with retained yeast fragments.
- MACS3 candidate regions are called at p <= 1e-3. Pooled matched-IgG peaks
  retained downstream require p <= 1e-4, q <= 0.01, and fold enrichment >= 3.
- Promoters are GENCODE M25 TSS +/- 1,000 bp. A peak must overlap a promoter by
  at least 200 bp. Final binding requires a promoter factor/IgG ratio >= 1.5 in
  both biological replicates.
- Whole-genome atomic peak segments retain their separate >=250-bp minimum.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.

To update CUT&RUN and regulatory outputs without rerunning RNA-seq or changing
bigWigs, use the existing candidate peak calls:

```bash
python 'pipeline/CUT&RUN/04_call_peaks_and_define_binding.py' --from-peaks
python 'pipeline/CUT&RUN/05_render_figures.py'
python pipeline/regulatory_target/06_define_regulatory_targets.py
python pipeline/regulatory_target/07_validate_outputs.py
```

Step 04 reapplies the configured criteria rather than reusing old promoter
memberships. Current promoter-bound totals are MCM3 2,945, NONO 745, and
PSPC1 4,903; regulatory-target totals are 962, 39, and 158, respectively.
Regulatory groups A/B/C/D contain 5/40/2/11 genes. These are threshold-dependent
results, not evidence that the chosen cutoffs are independently optimal.

The threshold update retained the existing coverage tracks. The historical
yeast scaling represented in these tracks differs from the current yeast-count
manifest; this separate normalization-provenance issue remains unresolved.
See `docs/Promoter_threshold_update.md` and the CUT&RUN source workbook.

### Binding Groups In Gene Profiles

TSS/TES profiles use the final promoter-bound gene sets, with each gene counted
once. The MCM3-containing groups are mutually exclusive: all three factors,
MCM3-NONO only, MCM3-PSPC1 only, and MCM3 only. "Only" describes which factors
pass the complete binding criteria; it does not establish absence of the other
factors. Promoter co-binding does not require identical peak coordinates or
demonstrate simultaneous occupancy in the same cells.

The retained filenames `MCM3_Peak_profiles.png` and
`MCM3_NONO_PSPC1_peak_profiles.png` are gene-level TSS/TES profiles, not
peak-centered measurements. Their gene lists are the corresponding
`cutrun_work/data/*_genes.tsv` files. Peak-promoter assignment tables are
interval annotations only; a gene can be assigned to multiple peak-segment
classes and must not inherit a gene-specific label from any single segment.
`Venn_Peaks.png` remains an interval-level analysis.

To regenerate only these two profiles and their source data:

```bash
python 'pipeline/CUT&RUN/05_render_figures.py' --promoter-groups-only
```

## Outputs

- `deg_work/`: RNA-seq quantifications, data tables, and figures. Salmon output
  is under `salmon`; native Salmon files such as `quant.sf` retain their
  standard names. Publication-compatible outputs
  are organized as `limma_outputs/DEGs`, `limma_outputs/heatmaps`, and
  `limma_outputs/volcanos`.
- `cutrun_work/`: alignments, BAMs, bigWigs, peaks, promoter sets, and figures.
  MACS peak calls are under `04_peaks`; final promoter evidence is under
  `05_promoters`.
  Final figures are in `visuals`; generated plotting inputs are in
  `data/figure_inputs`. The four publication-viewing tracks are directly under
  `03_bigwig/IGV_representation`.
- `regulatory_work/`: regulatory-target tables and figures.
- `docs/output_manifest.tsv`: complete final-figure inventory.
- `docs/data_manifest.tsv`: principal data-table inventory.

The numbered scripts and assay-specific instructions are in `pipeline/`.
