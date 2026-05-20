options(warn = -1)
suppressPackageStartupMessages(library('coloc'))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
type <- args[2]
tissue <- args[3]
locus_name <- args[4]
n_gwa <- as.integer(args[5])
n_qtl <- as.integer(args[6])
input_path <- args[7]
output_path <- args[8]
cols_save <- unlist(strsplit(args[9], ","))

result_cols <- c(c('Prefix','Locus','Tissue','Phenotype_id','Gene','Gene_id','TSS_distance',
					'SNP','REF','ALT','BETA_QTL','P_GWA','P_QTL','SNP.PP.H4',
					'PP.H0','PP.H1','PP.H2','PP.H3','PP.H4'),cols_save)
df_result <- as.data.frame(matrix(nrow=0,ncol=length(result_cols)))
names(df_result) <- result_cols

input <- read.table(file=sprintf('%s/%s_%s_%s.csv',input_path,prefix,locus_name,tissue), header=T,sep=',')
input <- na.omit(input)

pheno_lst <- unique(input$phenotype_id)
for (pheno in pheno_lst)
{
  temp <- input[which(input$phenotype_id == pheno),]
  gene <- temp[1,'gene_name']
  gene_id <- temp[1,'gene_id']
  if (type=='cc'){
  result <- coloc.abf(dataset1=list(snp=temp$SNP,pvalues=temp$P,type='cc',s=0.5,N=n_gwa), 
                      dataset2=list(snp=temp$SNP,pvalues=temp$pval_nominal,type='quant',N=n_qtl),MAF=temp$MAF)
				 }
  else if (type=='quant'){
  result <- coloc.abf(dataset1=list(snp=temp$SNP,pvalues=temp$P,type='quant',N=n_gwa), 
                      dataset2=list(snp=temp$SNP,pvalues=temp$pval_nominal,type='quant',N=n_qtl),MAF=temp$MAF)
                         }
  else {next}
  n <- result$summary['nsnps'][[1]]
  h0 <- result$summary['PP.H0.abf'][[1]]
  h1 <- result$summary['PP.H1.abf'][[1]]
  h2 <- result$summary['PP.H2.abf'][[1]]
  h3 <- result$summary['PP.H3.abf'][[1]]
  h4 <- result$summary['PP.H4.abf'][[1]]
  df_results <- result$results
  t <- df_results[order(df_results$SNP.PP.H4,decreasing = 'True'),][1,c('snp','SNP.PP.H4')]
  lead <- t$snp
  snp_h4 <- t$SNP.PP.H4
  p_gwa <- temp[which(temp$SNP == lead),'P']
  p_qtl <- temp[which(temp$SNP == lead),'pval_nominal']
  tss <- temp[which(temp$SNP == lead),'tss_distance']
  beta_qtl <- temp[which(temp$SNP == lead),'slope']
  ref <- temp[which(temp$SNP == lead),'ref_gwa']
  alt <- temp[which(temp$SNP == lead),'alt_gwa']
  #save results
  c_temp <- list(Prefix=prefix,Locus=locus_name,Tissue=tissue,REF=ref,ALT=alt,BETA_QTL=beta_qtl,
                 Phenotype_id=pheno,Gene=gene,Gene_id=gene_id,TSS_distance=tss,SNP=lead,P_GWA=p_gwa,P_QTL=p_qtl,SNP.PP.H4=snp_h4,
                 PP.H0=h0,PP.H1=h1,PP.H2=h2,PP.H3=h3,PP.H4=h4)
  #for other additional columns
  for (var_name in cols_save)
  {
    value <- temp[which(temp$SNP == lead),var_name]
    c_temp[[var_name]] <- value
  }
  df_result <- rbind(df_result, data.frame(c_temp))
}


if (nrow(df_result) >= 1) {
  write.csv(df_result, sprintf('%s/%s_%s_%s.csv', output_path, prefix, locus_name, tissue), row.names=FALSE)
} else {
  cat(sprintf("No rows to write for %s %s %s",prefix, locus_name, tissue))
}


