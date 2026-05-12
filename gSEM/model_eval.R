library(GenomicSEM)
library(Matrix)
library(stats)

root_path <- "/path/to/genomicsem/sumstats/"

load(paste0(root_path, "LDSCoutput_v2.RData"))

#---------------------------------------common factor model
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
EFA<-factanal(covmat = Ssmooth, factors = 1, rotation = "promax", nstart=5)

model <- "
F1 =~ NA*bladder + breast + cervical + colorectal + endometrial + kidney + lung + melanoma + oral + ovarian + prostate + thyroid
F1 ~~ 1*F1
"
fit <- usermodel(LDSCoutput, model = model, estimation = "DWLS", imp_cov = TRUE)
print(fit$results)
print(fit$modelfit)


#---------------------------------------two factor model
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
EFA<-factanal(covmat = Ssmooth, factors = 2, rotation = "promax", nstart=5)

model <- "
F1 =~ NA*bladder + breast + endometrial + kidney + ovarian + colorectal + melanoma + prostate + thyroid
F2 =~ NA*cervical + lung + oral
F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~~ F2
"
fit <- usermodel(LDSCoutput, model = model, std.lv = FALSE, estimation = "DWLS")
fit <- usermodel(LDSCoutput, model = model, estimation = "DWLS", imp_cov = TRUE)
print(fit$results)
print(fit$modelfit)


#---------------------------------------three factor model
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax", nstart=5)
EFA$loadings

model <- "
F1 =~ NA*cervical + lung + oral
F2 =~ NA*breast + colorectal + endometrial + kidney + melanoma + ovarian + prostate + thyroid
F3 =~ NA*bladder
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
F1 ~~ F2
F1 ~~ F3
F2 ~~ F3
bladder ~~ 0*bladder
"
fit <- usermodel(LDSCoutput, model = model, std.lv = FALSE, estimation = "DWLS")
fit <- usermodel(LDSCoutput, model = model, estimation = "DWLS", imp_cov = TRUE)
print(fit$results)
print(fit$modelfit)
