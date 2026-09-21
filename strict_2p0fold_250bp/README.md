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
  at least 250 bp. Final binding requires a promoter factor/IgG ratio >= 2 in
  both biological replicates.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.

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
