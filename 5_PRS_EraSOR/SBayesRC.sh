#!/bin/bash
#SBATCH -p cn
#SBATCH -N 1
#SBATCH -n 80
#SBATCH -J PRS
#SBATCH --array=1-14%1
#SBATCH -o PRS_%a.o
#SBATCH -e PRS_%a.e

set -euo pipefail

traits=(bladder breast cervical colorectal endometrial kidney lung melanoma oral ovarian prostate thyroid meta pan)

idx=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $idx -lt 0 || $idx -ge ${#traits[@]} ]]; then
  echo "Array index $SLURM_ARRAY_TASK_ID is out of range"
  exit 1
fi

trait="${traits[$idx]}"
echo "Running trait: ${trait}"

base="/path/to/PRS_EraSOR"
ma_file="${base}/${trait}.ma"
ld_folder="/path/to/SBayesRC/ukbEUR_Imputed/"
annot="/path/to/SBayesRC/annot_baseline2.2.txt"
out_folder="${base}/${trait}"
out_prefix="${out_folder}/${trait}"
threads=76

if [[ ! -f "${ma_file}" ]]; then
  echo "Input file not found: ${ma_file}"
  exit 2
fi
mkdir -p "$out_folder"

export OMP_NUM_THREADS="${threads}"

Rscript -e "SBayesRC::tidy(mafile='${ma_file}', LDdir='${ld_folder}', output='${out_prefix}_tidy.ma', log2file=TRUE)"
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='${ld_folder}', output='${out_prefix}_imp.ma', log2file=TRUE)"
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='${ld_folder}', outPrefix='${out_prefix}_sbrc', annot='${annot}', log2file=TRUE)"
