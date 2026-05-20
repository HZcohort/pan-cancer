#!/bin/bash
#SBATCH -p cn
#SBATCH -N 1
#SBATCH -n 80
#SBATCH -J partitioned_ldsc
#SBATCH -o logs/partitioned_ldsc.o
#SBATCH -e logs/partitioned_ldsc.e

source activate

LDSC_DIR="/path/to/ldsc"
PAN_SUMSTATS="/path/to/pan.gz"
HM3_SNPLIST="/path/to/w_hm3.snplist"
PARTITION_LDSC_DIR="/path/to/partitioned_ldscore"
OUT_DIR="/path/to/output/partitioned_ldsc"

DATA_DIR="${OUT_DIR}/data"
RESULT_DIR="${OUT_DIR}/result"
LOG_FILE="${OUT_DIR}/log.txt"

BASELINE_LD_V2="${PARTITION_LDSC_DIR}/EUR/baselineLD_v2/baselineLD."
BASELINE="${PARTITION_LDSC_DIR}/EUR/baseline/baseline."
WEIGHTS="${PARTITION_LDSC_DIR}/EUR/weights_hm3_no_hla/weights."
MULTI_TISSUE_EXP="${PARTITION_LDSC_DIR}/EUR/Multi_tissue_exp"
MULTI_TISSUE_CHROMATIN="${PARTITION_LDSC_DIR}/EUR/Multi_tissue_chromatin"

start_time=$(date +%s)
mkdir -p "${DATA_DIR}" "${RESULT_DIR}"
echo "Partitioned LDSC: start	($(date))" >> "${LOG_FILE}"

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

python "${LDSC_DIR}/ldsc.py" \
  --h2 "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE_LD_V2}" \
  --w-ld-chr "${WEIGHTS}" \
  --overlap-annot \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_baselineLD" &
wait
sleep 1

(cd "${MULTI_TISSUE_EXP}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts GTEx_tissue_exp.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_GTEx_tissue_exp") &

(cd "${MULTI_TISSUE_EXP}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts Frankelab_tissue_exp.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_Frankelab_tissue_exp") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_DNase.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_DNase") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_H3K4me1.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_H3K4me1") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_H3K4me3.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_H3K4me3") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_H3K9ac.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_H3K9ac") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_H3K27ac.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_H3K27ac") &

(cd "${MULTI_TISSUE_CHROMATIN}" && python "${LDSC_DIR}/ldsc.py" \
  --h2-cts "${DATA_DIR}/pan.sumstats.gz" \
  --ref-ld-chr "${BASELINE}" \
  --w-ld-chr "${WEIGHTS}" \
  --ref-ld-chr-cts chromatin_H3K36me3.ldcts \
  --not-M-5-50 \
  --out "${RESULT_DIR}/pan_chromatin_H3K36me3") &
wait
sleep 1

end_time=$(date +%s)
time_difference=$(((end_time - start_time) / 60))
echo "Partitioned LDSC: finished; Time elapsed: ${time_difference} minutes	($(date))" >> "${LOG_FILE}"
