SCHEME STDERR
MARKER snp
ALLELE A1 A2
EFFECT beta
STDERR se
PVALUE p

PROCESS /path/to/erasor_adjusted_gwas/bladder_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/breast_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/cervical_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/colorectal_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/endometrial_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/kidney_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/lung_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/melanoma_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/oral_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/ovarian_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/prostate_adjusted.txt
PROCESS /path/to/erasor_adjusted_gwas/thyroid_adjusted.txt

OUTFILE /path/to/PRS_EraSOR/meta_ .tbl
ANALYZE HETEROGENEITY
QUIT
