# MCM3 Project Analysis Pipeline

This directory contains the RNA-seq, CUT&RUN, and regulatory-target analyses
used for the final figures. All final figures are generated inside this
directory; no pipeline stage copies a figure from another reanalysis folder.

## Run The Workflow

Create the software environment:

```bash
mamba env create -f environment.yml
conda activate mcm3-final
```

Rebuild tables and figures from the local accepted alignment, peak, and Salmon
intermediates:

```bash
python pipeline/run_all.py --source accepted --step all --threads 16
```

Rebuild from raw FASTQs, including RNA-seq quantification and CUT&RUN
alignment, filtering, track generation, and peak calling:

```bash
python pipeline/run_all.py --source raw --step all --threads 16
```

Raw FASTQs used by this workflow are stored beside this directory in
`bulkRNAseq_data/` and `CUT&RUN_data/`. CUT&RUN files are organized under
`20191012/`, `20200923/`, and `20200929/`; archived cross-knockdown and
undetermined reads are excluded from analysis by the generated sample manifest.

Use `--dry-run` to print the stages without running them.

## Analysis Summary

- RNA-seq reads are quantified with Salmon and tested with limma-voom. DEGs
  require BH-adjusted p-value <= 0.05 and absolute log2 fold change >= 0.28.
- CUT&RUN retains properly paired, primary, non-duplicate mouse alignments with
  MAPQ >= 30. Coverage tracks are normalized with retained yeast fragments.
- Peaks require p <= 1e-4 and fold enrichment >= 3. The q-value cutoff is
  <= 0.05 for 2019 and <= 0.01 for 2020.
- Promoters are GENCODE M25 TSS +/- 1,000 bp. A peak must overlap a promoter by
  at least 250 bp, and final binding requires support in at least 4 of 6 factor
  libraries. The 2020 factor/IgG ratio must be >= 2 in both replicates.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.

## Outputs

- `deg_work/`: RNA-seq quantifications, data tables, and figures. Salmon output
  is under `salmon`; native Salmon files such as `quant.sf` retain their
  standard names. Publication-compatible outputs
  are organized as `limma_outputs/DEGs`, `limma_outputs/heatmaps`, and
  `limma_outputs/volcanos`.
- `cutrun_work/`: alignments, BAMs, bigWigs, peaks, promoter sets, and figures.
  MACS peak calls from both years are unified under `04_peaks/2019` and
  `04_peaks/2020`; final promoter evidence is under `05_promoters`.
  Final figures are in `visuals`; generated plotting inputs are in
  `data/figure_inputs`. The four publication-viewing tracks are directly under
  `03_bigwig/IGV_representation`.
- `regulatory_work/`: regulatory-target tables and figures.
- `docs/output_manifest.tsv`: complete final-figure inventory.
- `docs/data_manifest.tsv`: principal data-table inventory.

The numbered scripts and assay-specific instructions are in `pipeline/`.
