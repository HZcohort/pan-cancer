#!/bin/bash
#SBATCH -p cn
#SBATCH -N 1
#SBATCH -n 80
#SBATCH -J phewas_gc
#SBATCH -o logs/phewas_gc.o
#SBATCH -e logs/phewas_gc.e

source activate

LDSC_DIR="/path/to/ldsc"
PAN_SUMSTATS="/path/to/pan.gz"
HM3_SNPLIST="/path/to/w_hm3.snplist"
TRAIT_SUMSTATS_DIR="/path/to/phenomewide_EUR"
LDSC_REF_DIR="/path/to/ld_score/EUR"
OUT_DIR="/path/to/output/phewas_gc"
MERGE_SCRIPT="/path/to/ldsc_merge.py"

DATA_DIR="${OUT_DIR}/data"
RESULT_TRAITS_DIR="${OUT_DIR}/result_traits"
LOG_FILE="${OUT_DIR}/log.txt"

start_time=$(date +%s)
mkdir -p "${DATA_DIR}" "${RESULT_TRAITS_DIR}"
echo "Phenome-wide genetic correlation estimation: start	($(date))" >> "${LOG_FILE}"

conda activate ldsc
python "${LDSC_DIR}/munge_sumstats.py" \
  --sumstats "${PAN_SUMSTATS}" \
  --merge-alleles "${HM3_SNPLIST}" \
  --out "${DATA_DIR}/pan" \
  --snp snp \
  --p p \
  --a1 A1 \
  --a2 A2 \
  --N-col N \
  --signed-sumstats beta,0 \
  --chunksize 500000 &
wait
sleep 1

conda activate ldsc
cd "${TRAIT_SUMSTATS_DIR}"
files=$(ls *.sumstats.gz | tr "\n" "," | sed "s/,$//")
python "${LDSC_DIR}/ldsc.py" \
  --rg "${DATA_DIR}/pan.sumstats.gz,${files}" \
  --ref-ld-chr "${LDSC_REF_DIR}/" \
  --w-ld-chr "${LDSC_REF_DIR}/" \
  --out "${RESULT_TRAITS_DIR}/ldsc_pan_traits" &
wait
sleep 1

conda activate base
python "${MERGE_SCRIPT}" \
  --input_path "${RESULT_TRAITS_DIR}" \
  --output_path "${OUT_DIR}" \
  --prefix traits &
wait
sleep 1

end_time=$(date +%s)
time_difference=$(((end_time - start_time) / 60))
echo "Phenome-wide genetic correlation estimation: finished; Time elapsed: ${time_difference} minutes	($(date))" >> "${LOG_FILE}"
