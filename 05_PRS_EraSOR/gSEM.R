library(GenomicSEM)
library(Matrix)
library(stats)

root_path <- "/path/to/erasor_adjusted_gwas/"

traits <- c(paste0(root_path, "bladder_adjusted.gz"),
            paste0(root_path, "breast_adjusted.gz"),
            paste0(root_path, "cervical_adjusted.gz"),
            paste0(root_path, "colorectal_adjusted.gz"),
            paste0(root_path, "endometrial_adjusted.gz"),
            paste0(root_path, "kidney_adjusted.gz"),
            paste0(root_path, "lung_adjusted.gz"),
            paste0(root_path, "melanoma_adjusted.gz"),
            paste0(root_path, "oral_adjusted.gz"),
            paste0(root_path, "ovarian_adjusted.gz"),
            paste0(root_path, "prostate_adjusted.gz"),
            paste0(root_path, "thyroid_adjusted.gz"))

trait.names <- c("bladder", "breast", "cervical", "colorectal", "endometrial",
                 "kidney", "lung", "melanoma", "oral", "ovarian", "prostate", "thyroid")

sample.prev <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

population.prev <- c(
  0.021,
  0.130,
  0.006,
  0.039,
  0.031,
  0.018,
  0.052,
  0.022,
  0.012,
  0.011,
  0.132,
  0.011
)

ld <- "/path/to/ld_score/LDSC/EUR/"
wld <- "/path/to/ld_score/LDSC/EUR/"

LDSCoutput <- ldsc(traits = traits,
                   sample.prev = sample.prev,
                   population.prev = population.prev,
                   ld = ld,
                   wld = wld,
                   trait.names = trait.names)

save(LDSCoutput, file = paste0(root_path, "LDSCoutput.RData"))

Ssmooth <- as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
EFA <- factanal(covmat = Ssmooth, factors = 3, rotation = "promax", nstart = 5)
EFA$loadings

N <- c(66001, 430306, 46708, 341185, 40250, 119798,
       192344, 190241, 48240, 104336, 383422, 32070)
se.logit <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
              TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

summary_path <- "/path/to/erasor_adjusted_gwas/"
files <- c(paste0(summary_path, "bladder_adjusted.txt"),
           paste0(summary_path, "breast_adjusted.txt"),
           paste0(summary_path, "cervical_adjusted.txt"),
           paste0(summary_path, "colorectal_adjusted.txt"),
           paste0(summary_path, "endometrial_adjusted.txt"),
           paste0(summary_path, "kidney_adjusted.txt"),
           paste0(summary_path, "lung_adjusted.txt"),
           paste0(summary_path, "melanoma_adjusted.txt"),
           paste0(summary_path, "oral_adjusted.txt"),
           paste0(summary_path, "ovarian_adjusted.txt"),
           paste0(summary_path, "prostate_adjusted.txt"),
           paste0(summary_path, "thyroid_adjusted.txt"))

ref <- "/path/to/reference.1000G.maf.0.005.txt"
Hail <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
          FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

cancer_sumstats <- sumstats(files = files,
                            ref = ref,
                            trait.names = trait.names,
                            se.logit = se.logit,
                            OLS = Hail,
                            linprob = Hail,
                            N = N,
                            parallel = TRUE,
                            cores = 14)

save(cancer_sumstats, file = paste0(root_path, "cancer_sumstats.RData"))

model <- "F1 =~ NA*bladder + breast + cervical + colorectal + endometrial + kidney + lung + melanoma + oral + ovarian + prostate + thyroid
F1~~1*F1"

fit <- usermodel(LDSCoutput,
                 estimation = "DWLS",
                 model = model,
                 CFIcalc = TRUE,
                 std.lv = TRUE,
                 imp_cov = FALSE)
save(fit, file = paste0(root_path, "one_factor_fit.RData"))

model <- "F1 =~ NA*bladder + breast + cervical + colorectal + endometrial + kidney + lung + melanoma + oral + ovarian + prostate + thyroid
F1~~1*F1
F1 ~ SNP"

pan_gwas <- userGWAS(covstruc = LDSCoutput,
                     SNPs = cancer_sumstats,
                     model = model,
                     sub = c("F1~SNP"),
                     smooth_check = TRUE,
                     fix_measurement = TRUE,
                     Q_SNP = TRUE,
                     cores = 10)

save(pan_gwas, file = paste0(root_path, "pan_eras_userGWAS.RData"))
write.table(pan_gwas, file = paste0(root_path, "pan_eras.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
