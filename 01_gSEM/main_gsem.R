library(GenomicSEM)
library(Matrix)
library(stats)

# -------------------------------------------------------------------------
# User-defined paths
# -------------------------------------------------------------------------
root_path <- "/path/to/genomicsem/sumstats/"
summary_path <- "/path/to/processed/gwas/"
ld <- "/path/to/LDSC/EUR/"
wld <- "/path/to/LDSC/EUR/"
ref <- "/path/to/reference.1000G.maf.0.005.txt"

# -------------------------------------------------------------------------
# Step 1: estimate LDSC genetic covariance and sampling covariance matrices
# -------------------------------------------------------------------------
traits <- c(paste0(root_path, "bladder.sumstats.gz"),
           paste0(root_path, "breast.sumstats.gz"),
           paste0(root_path, "cervical.sumstats.gz"),
           paste0(root_path, "colorectal.sumstats.gz"),
           paste0(root_path, "endometrial.sumstats.gz"),
           paste0(root_path, "kidney.sumstats.gz"),
           paste0(root_path, "lung.sumstats.gz"),
           paste0(root_path, "melanoma.sumstats.gz"),
           paste0(root_path, "oral.sumstats.gz"),
           paste0(root_path, "ovarian.sumstats.gz"),
           paste0(root_path, "prostate.sumstats.gz"),
           paste0(root_path, "thyroid.sumstats.gz"))

trait.names <- c("bladder", "breast", "cervical", "colorectal", "endometrial",
"kidney", "lung", "melanoma", "oral", "ovarian", "prostate", "thyroid")

sample.prev <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

population.prev <- c(
    0.021, # Bladder
    0.130, # Breast
    0.006, # Cervical
    0.039, # Colorectal
    0.031, # Endometrial
    0.018, # Kidney
    0.052, # Lung
    0.022, # Melanoma
    0.012, # Oral cavity and pharynx
    0.011, # Ovarian
    0.132, # Prostate
    0.011  # Thyroid
  )

LDSCoutput <- ldsc(traits = traits,
                   sample.prev = sample.prev,
                   population.prev = population.prev,
                   ld = ld,
                   wld = wld,
                   trait.names = trait.names)

save(LDSCoutput, file = paste0(root_path, "LDSCoutput_v2.RData"))

# -------------------------------------------------------------------------
# Step 2: process SNP-level summary statistics for GenomicSEM
# -------------------------------------------------------------------------
N <- c(69049, 430306, 47740, 354068, 46158, 123750, 200824, 204496, 48864, 108144, 406587, 33586)
se.logit <- c(T,T,T,T,T,T,T,T,T,T,T,T)

files <- c(paste0(summary_path, "bladder.txt"),
           paste0(summary_path, "breast.txt"),
           paste0(summary_path, "cervical.txt"),
           paste0(summary_path, "colorectal.txt"),
           paste0(summary_path, "endometrial.txt"),
           paste0(summary_path, "kidney.txt"),
           paste0(summary_path, "lung.txt"),
           paste0(summary_path, "melanoma.txt"),
           paste0(summary_path, "oral.txt"),
           paste0(summary_path, "ovarian.txt"),
           paste0(summary_path, "prostate.txt"),
           paste0(summary_path, "thyroid.txt"))

trait.names <- c("bladder", "breast", "cervical", "colorectal", "endometrial",
"kidney", "lung", "melanoma", "oral", "ovarian", "prostate", "thyroid")

Hail <- c(F,F,F,F,F,F,F,F,F,F,F,F)

cancer_sumstats <- sumstats(files=files, ref=ref, trait.names=trait.names,
se.logit=se.logit, OLS=Hail, linprob=Hail, N=N, parallel=TRUE, cores=14)

save(cancer_sumstats, file = paste0(root_path, "cancer_sumstats_v2.RData"))

# -------------------------------------------------------------------------
# Step 3: chromosome-wise multivariate GWAS of the pan-cancer factor
# -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  load(paste0(root_path, "LDSCoutput_v2.RData"))
  load(paste0(root_path, "cancer_sumstats_v2.RData"))

  chr_input <- as.numeric(args[1])
  if (is.na(chr_input)) {
    stop("Provided chromosome number is not valid. Please enter a numeric value.")
  }

  cancer_sumstats <- subset(cancer_sumstats, CHR == chr_input)
  cat("Running analysis for chromosome:", chr_input, "\n")

  model<-"F1 =~ NA*bladder + breast + cervical + colorectal + endometrial + kidney + lung + melanoma + oral + ovarian + prostate + thyroid
F1~~1*F1
F1 ~ SNP"

  CorrelatedFactors <- userGWAS(covstruc = LDSCoutput,
                                SNPs = cancer_sumstats,
                                model = model,
                                sub = c("F1~SNP"),
                                SNPSE=0.005,
                                smooth_check = TRUE,
                                fix_measurement = TRUE,
                                Q_SNP = TRUE,
                                cores = 79)

  save(CorrelatedFactors, file = paste0(root_path, "cancer_results_chr", chr_input, ".RData"))
}
