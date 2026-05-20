library(HDL)
library(GenomicSEM)
library(Matrix)
library(stats)
library(readr)
library(stringr)
library(lavaan)
library(dplyr)

root_path <- "/path/to/genomicsem/sumstats/"
hdl_ld_path <- "/path/to/HDL/UKB_imputed_hapmap2_SVD_eigen99_extraction/"
hdl_wrapper <- "/path/to/hdl.R"

if (file.exists(hdl_wrapper)) {
  source(hdl_wrapper)
}

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

trait.names <- c("bladder", "breast", "cervical", "colorectal", "endometrial",
"kidney", "lung", "melanoma", "oral", "ovarian", "prostate", "thyroid")

hdl.covstruct <- hdl(traits, sample.prev = sample.prev, population.prev = population.prev,
trait.names=trait.names, LD.path=hdl_ld_path, method = "piecewise")
saveRDS(hdl.covstruct, file=paste0(root_path, "hdl_rg.rds"))

hdl.covstruct <- readRDS(paste0(root_path, "hdl_rg.rds"))
model <- "
F1 =~ NA*bladder + breast + cervical + colorectal + endometrial + kidney + lung + melanoma + oral + ovarian + prostate + thyroid
F1 ~~ 1*F1
"
fit <- usermodel(hdl.covstruct, model = model, estimation = "DWLS", imp_cov = TRUE)
print(fit$results)
print(fit$modelfit)
