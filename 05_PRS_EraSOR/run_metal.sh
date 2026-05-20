#!/bin/bash
#SBATCH -p cn
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J metal_meta_prs
#SBATCH -o logs/metal_meta_prs.o
#SBATCH -e logs/metal_meta_prs.e

set -euo pipefail

METAL="/path/to/metal"
METAL_SCRIPT="/path/to/PRS_EraSOR/metal_meta_analysis.metal"

"${METAL}" "${METAL_SCRIPT}"
