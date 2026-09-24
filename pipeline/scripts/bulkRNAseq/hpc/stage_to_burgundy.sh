#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
REMOTE_HOST="${MCM3_HPC_HOST:-ysuzuki2@burgundy.hpc.cityu.edu.hk}"
REMOTE_ROOT="${MCM3_HPC_PROJECT_ROOT:-/scratch/ysuzuki2/MCM3_project_analysis_pipeline}"
SSH_COMMAND="ssh -o ProxyJump=none -i ${HOME}/.ssh/id_ed25519_axis4_burgundy -o IdentitiesOnly=yes"
FILE_LIST="$(mktemp)"
trap 'rm -f "${FILE_LIST}"' EXIT

python3 - "${PROJECT_ROOT}" > "${FILE_LIST}" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1]).resolve()
with (root / "pipeline/deg_work/metadata/Samples.tsv").open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["selected"].strip().lower() == "true":
            for column in ("r1", "r2"):
                print(pathlib.Path(row[column]).resolve().relative_to(root))
PY

${SSH_COMMAND} "${REMOTE_HOST}" "mkdir -p '${REMOTE_ROOT}/pipeline/reference' '${REMOTE_ROOT}/bulkRNAseq_data' '${REMOTE_ROOT}/logs'"

rsync -a --progress \
  --exclude '._*' --exclude '.cache' --exclude 'reference' --exclude 'cutrun_work' \
  --exclude 'regulatory_work' --exclude 'deg_work/igv_tracks' \
  -e "${SSH_COMMAND}" \
  "${PROJECT_ROOT}/pipeline/" "${REMOTE_HOST}:${REMOTE_ROOT}/pipeline/"

rsync -a --progress -e "${SSH_COMMAND}" \
  "${PROJECT_ROOT}/pipeline/reference/GRCm38.primary_assembly.gencodeM25_contigs.fa.gz" \
  "${PROJECT_ROOT}/pipeline/reference/gencode.vM25.annotation.gtf.gz" \
  "${PROJECT_ROOT}/pipeline/reference/tx2gene.tsv" \
  "${REMOTE_HOST}:${REMOTE_ROOT}/pipeline/reference/"

rsync -a --partial --progress --files-from="${FILE_LIST}" -e "${SSH_COMMAND}" \
  "${PROJECT_ROOT}/" "${REMOTE_HOST}:${REMOTE_ROOT}/"
