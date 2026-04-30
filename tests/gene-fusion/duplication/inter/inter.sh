#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${SCRIPT_DIR}/run_${RUN_TS}"
mkdir -p "${RUN_DIR}"
LOG="${RUN_DIR}/log.txt"

cd "${REPO_ROOT}"

# fixed inputs
IGENOME="examples/genomes/carsonella_ruddii_JRPAMB4.peg"
PSUB="psubmatrix.txt"
NGENOMES=3
RSEED=42

# fixed evolution params (everything off except fusion/sub-dup/intra)
GENE_VAR=0
LOC_VAR=0
GENE_DUP=0
GSET_VAR=0
GENE_ADD=0
GENE_FUSION=1.0
SUB_GENE_DUP=1.0
INTRA_SUB_GENE_DUP=0.0   # 0.0 -> "same gene" branch (line 799); 1.0 -> "different gene" branch (line 850)

# per-run files
HGT_POOL="${RUN_DIR}/hgt.pool"
TREE="${RUN_DIR}/run.tree"
OPREFIX="${RUN_DIR}/run"

{
echo "================================================================================"
echo "PANPROVA gene-fusion / sub-gene-duplication / inter"
echo "================================================================================"
echo "repo:   ${REPO_ROOT}"
echo "outdir: ${RUN_DIR}"
echo

# empty HGT pool (fusion test does not need HGT)
: > "${HGT_POOL}"

echo "================================================================================"
echo "Generating tree (${NGENOMES} nodes, seed ${RSEED})..."
cmd="python3 generate_tree_2.py ${NGENOMES} ${RSEED} ${TREE}"
echo "$cmd"
$cmd

echo
echo "================================================================================"
echo "Evolving (calling ./evolve directly, skipping PANPROVA.sh post-processing)..."
cmd="./evolve \
${IGENOME} \
${HGT_POOL} \
${OPREFIX} \
${TREE} \
${PSUB} \
${GENE_VAR} ${LOC_VAR} ${GENE_DUP} \
${GSET_VAR} ${GENE_ADD} \
${RSEED} \
${GENE_FUSION} ${SUB_GENE_DUP} ${INTRA_SUB_GENE_DUP}"
echo "$cmd"
date
$cmd
date

echo
echo "================================================================================"
echo "output:        ${RUN_DIR}"
echo "chimeras csv:  ${OPREFIX}.chimeras.csv"
echo "chimeras tsv:  ${OPREFIX}.chimeras.tsv"
} 2>&1 | tee "${LOG}"
