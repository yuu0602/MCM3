#!/usr/bin/env bash
#SBATCH --job-name=mcm3-rnaseq-igv
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=2-12:00:00
#SBATCH --output=/scratch/%u/MCM3_project_analysis_pipeline/logs/igv_tracks_%j.out
#SBATCH --error=/scratch/%u/MCM3_project_analysis_pipeline/logs/igv_tracks_%j.err

set -euo pipefail

PROJECT_ROOT="${MCM3_HPC_PROJECT_ROOT:-/scratch/${USER}/MCM3_project_analysis_pipeline}"
ENVIRONMENT="${MCM3_HPC_ENVIRONMENT:-/scratch/${USER}/mcm3-pipeline-env}"
MICROMAMBA="${MCM3_HPC_MICROMAMBA:-/scratch/${USER}/micromamba/bin/micromamba}"
NODE_ROOT="${SLURM_TMPDIR:-/tmp/${USER}/mcm3_rnaseq_igv_${SLURM_JOB_ID}}"
NODE_PROJECT="${NODE_ROOT}/MCM3_project_analysis_pipeline"

mkdir -p "${NODE_ROOT}"
trap 'rm -rf "${NODE_ROOT}"' EXIT

echo "[$(date)] Staging RNA-seq pipeline and selected FASTQs to ${NODE_ROOT}"
rsync -a --delete --exclude 'cutrun_work' --exclude 'deg_work/igv_tracks' "${PROJECT_ROOT}/pipeline/" "${NODE_PROJECT}/pipeline/"
rsync -a "${PROJECT_ROOT}/bulkRNAseq_data/" "${NODE_PROJECT}/bulkRNAseq_data/"

gzip -dc "${NODE_PROJECT}/pipeline/reference/GRCm38.primary_assembly.gencodeM25_contigs.fa.gz" > "${NODE_PROJECT}/pipeline/reference/GRCm38.primary_assembly.gencodeM25_contigs.fa"
gzip -dc "${NODE_PROJECT}/pipeline/reference/gencode.vM25.annotation.gtf.gz" > "${NODE_PROJECT}/pipeline/reference/gencode.vM25.annotation.gtf"

sed -i "s|${PROJECT_ROOT}|${NODE_PROJECT}|g" "${NODE_PROJECT}/pipeline/deg_work/metadata/Samples.tsv"

export MCM3_HISAT2_INDEX_TMPDIR="${NODE_ROOT}/hisat2_index"
export MPLCONFIGDIR="${NODE_ROOT}/matplotlib"

echo "[$(date)] Generating RNA-seq IGV tracks"
"${MICROMAMBA}" run -p "${ENVIRONMENT}" python "${NODE_PROJECT}/pipeline/scripts/bulkRNAseq/03_render_figures.py" \
  --igv-tracks \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] Returning bigWigs and manifest to shared storage"
mkdir -p "${PROJECT_ROOT}/pipeline/deg_work/igv_tracks"
rsync -a --delete --exclude 'alignments' \
  "${NODE_PROJECT}/pipeline/deg_work/igv_tracks/" \
  "${PROJECT_ROOT}/pipeline/deg_work/igv_tracks/"

echo "[$(date)] Complete"
