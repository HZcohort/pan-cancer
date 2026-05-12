#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/path/to/pan_cancer_gwas"
SCRIPT_DIR="/path/to/consensus_mapping_scripts"
RESULT_DIR="${BASE_DIR}/result"
GWAS_DIR="${BASE_DIR}/GWAS/processed"
ANNOTATION_DIR="/path/to/annotation"
REF_PREFIX="/path/to/1000G/EUR/g1000_EUR"
POPS_DIR="/path/to/pops"
N_GWA="6127"
COLS_SAVE="ref,alt"

mkdir -p "${RESULT_DIR}/eQTL/result" "${RESULT_DIR}/sQTL/result" "${RESULT_DIR}/pQTL/result"
mkdir -p "${RESULT_DIR}/PAV" "${RESULT_DIR}/H-MAGMA/data" "${RESULT_DIR}/PoPS/MAGMA" "${RESULT_DIR}/PoPS/data"
mkdir -p "${RESULT_DIR}/consensus"

# QTL colocalization assumes that locus-level harmonized GWAS-QTL input files have
# already been prepared in the corresponding data directories.
line_counter=0
while IFS=" " read -r prefix tissue loci n_gwa n_qtl type
do
    IFS="," read -ra loci_array <<< "${loci}"
    for locus in "${loci_array[@]}"
    do
        Rscript "${SCRIPT_DIR}/coloc_eqtl.R" "${prefix}" "${type}" "${tissue}" "${locus}" "${n_gwa}" "${n_qtl}" \
            "${RESULT_DIR}/eQTL/data" "${RESULT_DIR}/eQTL/result" "${COLS_SAVE}" &
        ((line_counter+=1))
        if ((line_counter % 30 == 0)); then wait; fi
    done
done < "${RESULT_DIR}/eQTL/eQTL.parameter"
wait
find "${RESULT_DIR}/eQTL/result" -name '*.csv' -exec awk 'NR == 1 || FNR > 1' {} + > "${RESULT_DIR}/eQTL.csv"

line_counter=0
while IFS=" " read -r prefix tissue loci n_gwa n_qtl type
do
    IFS="," read -ra loci_array <<< "${loci}"
    for locus in "${loci_array[@]}"
    do
        Rscript "${SCRIPT_DIR}/coloc_sqtl.R" "${prefix}" "${type}" "${tissue}" "${locus}" "${n_gwa}" "${n_qtl}" \
            "${RESULT_DIR}/sQTL/data" "${RESULT_DIR}/sQTL/result" "${COLS_SAVE}" &
        ((line_counter+=1))
        if ((line_counter % 40 == 0)); then wait; fi
    done
done < "${RESULT_DIR}/sQTL/sQTL.parameter"
wait
find "${RESULT_DIR}/sQTL/result" -name '*.csv' -exec awk 'NR == 1 || FNR > 1' {} + > "${RESULT_DIR}/sQTL.csv"

line_counter=0
while IFS=" " read -r prefix tissue locus script data_path n_gwa n_qtl type
do
    Rscript "${SCRIPT_DIR}/coloc_pqtl.R" "${prefix}" "${type}" "${tissue}" "${locus}" "${n_gwa}" "${n_qtl}" \
        "${RESULT_DIR}/pQTL/data" "${RESULT_DIR}/pQTL/result" "${COLS_SAVE}" &
    ((line_counter+=1))
    if ((line_counter % 40 == 0)); then wait; fi
done < "${RESULT_DIR}/pQTL/pQTL.parameter"
wait
find "${RESULT_DIR}/pQTL/result" -name '*.csv' -exec awk 'NR == 1 || FNR > 1' {} + > "${RESULT_DIR}/pQTL.csv"

python "${SCRIPT_DIR}/pav.py" \
    --input_path "${RESULT_DIR}/loci_verified" \
    --output_path "${RESULT_DIR}/PAV" \
    --bim_file "${REF_PREFIX}.bim" \
    --flip_kb 2000

line_counter=0
while IFS=" " read -r prefix locus
do
    plink --bfile "${REF_PREFIX}" \
        --clump "${RESULT_DIR}/PAV/${prefix}_${locus}.indsig" "${RESULT_DIR}/PAV/${prefix}_${locus}.sumstats" \
        --clump-kb 2000 \
        --clump-r2 0.80 \
        --clump-allow-overlap \
        --clump-index-first \
        --out "${RESULT_DIR}/PAV/${prefix}_${locus}" &
    ((line_counter+=1))
    if ((line_counter % 40 == 0)); then wait; fi
done < "${RESULT_DIR}/PAV/pav.parameter"
wait

python "${SCRIPT_DIR}/clump_merge.py" --input_path "${RESULT_DIR}/PAV" --output_path "${RESULT_DIR}/PAV"
awk -F'\t' 'NR==FNR{filter[$1]; next} $1 in filter' "${RESULT_DIR}/PAV/proxy_snp.list" \
    "${ANNOTATION_DIR}/VEP/EUR_vep.tsv" > "${RESULT_DIR}/PAV/VEP_selected.tsv"
python "${SCRIPT_DIR}/pav_select.py" --input_path "${RESULT_DIR}/PAV" --output_path "${RESULT_DIR}"

python "${SCRIPT_DIR}/clinvar.py" \
    --input_path "${RESULT_DIR}/loci_verified" \
    --output_path "${RESULT_DIR}" \
    --clinvar_data_path "${ANNOTATION_DIR}/ClinVar" \
    --gene_annotation_file "${ANNOTATION_DIR}/gene/gene_annotation.txt" \
    --hpo_dict "{'pan': 'HP:0034930,HP:0004377,HP:0100001,HP:0200013,HP:0012288,HP:0009728,HP:0100013,HP:0012780,HP:0100568,HP:0100012,HP:0007379,HP:0100544,HP:0100604,HP:0100826,HP:0004375,HP:0100606,HP:0010622,HP:0008069,HP:0100521,HP:0031459,HP:0100742,HP:0002898,HP:0031492,HP:0012316,HP:0010566,HP:0002861,HP:0030060,HP:0034557,HP:0100242'}"

python "${SCRIPT_DIR}/mgi.py" \
    --input_path "${RESULT_DIR}/loci_verified" \
    --output_path "${RESULT_DIR}" \
    --mgi_data_path "${ANNOTATION_DIR}/MGI" \
    --gene_annotation_file "${ANNOTATION_DIR}/gene/gene_annotation.txt" \
    --mpo_dict "{'pan': 'MP:0000858,MP:0003448,MP:0021056,MP:0010537,MP:0002019,MP:0010307,MP:0002009'}"

zcat "${GWAS_DIR}/pan.gz" | awk 'BEGIN {FS=OFS="\t"} NR==1 {for(i=1;i<=NF;i++) header[$i]=i; print "snp", "p"} NR>1 {print $header["snp"], $header["p"]}' > "${RESULT_DIR}/H-MAGMA/data/pan"
line_counter=0
while IFS=" " read -r prefix tissue loci_lst n_gwa
do
    magma --bfile "${REF_PREFIX}" \
        --pval "${RESULT_DIR}/H-MAGMA/data/${prefix}" use=snp,p N="${n_gwa}" \
        --gene-annot "${ANNOTATION_DIR}/HMAGMA/${tissue}.transcript.annot" \
        --gene-model snp-wise=mean \
        --out "${RESULT_DIR}/H-MAGMA/data/${prefix}_${tissue}" &
    ((line_counter+=1))
    if ((line_counter % 40 == 0)); then wait; fi
done < "${RESULT_DIR}/H-MAGMA/HMAGMA.parameter"
wait

python "${SCRIPT_DIR}/HMAGMA_merge.py" \
    --parameter_file "${RESULT_DIR}/H-MAGMA/HMAGMA.parameter" \
    --gene_annot_file "${ANNOTATION_DIR}/HMAGMA/ENSG_symbol.csv" \
    --input_path "${RESULT_DIR}/H-MAGMA/data" \
    --output_path "${RESULT_DIR}"

awk 'BEGIN{FS=" "} {print $1,$3,$4}' "${RESULT_DIR}/H-MAGMA/HMAGMA.parameter" | sort | uniq > "${RESULT_DIR}/PoPS/MAGMA/PoPS.parameter"
while IFS=" " read -r prefix loci n_gwa
do
    magma --bfile "${REF_PREFIX}" \
        --pval "${RESULT_DIR}/H-MAGMA/data/${prefix}" use=snp,p N="${n_gwa}" \
        --gene-annot "${POPS_DIR}/ENSEMBLE.genes.annot" \
        --gene-model snp-wise=mean \
        --out "${RESULT_DIR}/PoPS/MAGMA/${prefix}"
done < "${RESULT_DIR}/PoPS/MAGMA/PoPS.parameter"

while IFS=" " read -r prefix loci n_gwa
do
    python "${POPS_DIR}/pops.py" \
        --gene_annot_path "${POPS_DIR}/gene_annot_jun10.txt" \
        --feature_mat_prefix "${POPS_DIR}/features_munged/pops_features" \
        --num_feature_chunks 123 \
        --magma_prefix "${RESULT_DIR}/PoPS/MAGMA/${prefix}" \
        --out_prefix "${RESULT_DIR}/PoPS/data/${prefix}" \
        --verbose
done < "${RESULT_DIR}/PoPS/MAGMA/PoPS.parameter"

python "${SCRIPT_DIR}/pops_merge.py" \
    --input_path "${RESULT_DIR}/PoPS" \
    --output_path "${RESULT_DIR}" \
    --gene_annot_file "${POPS_DIR}/gene_annot_jun10.txt"

CONSENSUS_DICT="$(tr -d '\n' < "${SCRIPT_DIR}/consensus_dict_pan.txt")"
python "${SCRIPT_DIR}/consensus.py" \
    --mode_GWAS \
    --input_path "${RESULT_DIR}" \
    --input_path_loci "${RESULT_DIR}/loci_verified" \
    --output_path "${RESULT_DIR}/consensus" \
    --output_path_final "${RESULT_DIR}" \
    --hmagma_name_file "${ANNOTATION_DIR}/HMAGMA/Cell_IDs.xlsx" \
    --gene_update_file "${ANNOTATION_DIR}/gene/update_old_new.tsv" \
    --gene_annotation_file "${ANNOTATION_DIR}/gene/ENSGID_NAME.tsv" \
    --consensus_dict "${CONSENSUS_DICT}" \
    --coloc_pp4 "[0.8, 0.8, 0.8]"
