# Pan-cancer Genetics Analysis Code

This repository contains analysis scripts for the manuscript "Genomic structural equation modelling of 12 cancers identifies a latent pan-cancer susceptibility factor and shared determinants of carcinogenesis". Each subfolder corresponds to one analysis component and includes a section-specific `README.md` describing the expected inputs, script order, and outputs.

## Manuscript Information

**Title:** Genomic structural equation modelling of 12 cancers identifies a latent pan-cancer susceptibility factor and shared determinants of carcinogenesis

**Authors:** Qian Xu, Lu Niu, Zian Cao, Jiwen Geng, Huazhen Yang, Yu Zeng, Yanan Zhang, Qian Li, Chen Suo, Gang Yuan, Huan Song, Can Hou

## Folder Overview

- `gSEM/`: GenomicSEM covariance estimation, model evaluation, LDSC and HDL comparison, and multivariate GWAS of the latent pan-cancer factor.
- `risk_loci_annotation/`: GWAS summary-statistic processing, LD-independent risk locus definition, lead variant verification, and GWAS Catalog novelty lookup.
- `partitioned_heritability/`: Stratified LDSC analyses for baseline, tissue, cell-type, histone-mark, and expression annotations.
- `phenowide_gc/`: LDSC-based genetic correlations between the pan-cancer factor and external GWAS traits.
- `PRS_EraSOR/`: EraSOR adjustment, constructed pan-cancer GWAS, cross-cancer meta-analysis PRS, SBayesRC PRS weights, P+T PRS generation, and PLINK scoring preparation.
- `PRS_evaluation/`: UK Biobank prospective PRS evaluation and East Asian case-control validation.
- `phewas/`: UK Biobank PheWAS of the pan-cancer PRS using DiNetxify.
- `consensus_mapping/`: Consensus-based gene mapping using QTL colocalization, PAV linkage, ClinVar, MGI, H-MAGMA, PoPS, and final evidence aggregation.
- `pathway_enrichment/`: MAGMA gene-set enrichment and mapped-gene enrichment analyses.
- `pQTL_MR/`: Proteome-wide cis-pQTL MR, replicated protein selection, and SMR/HEIDI verification workflow.

## Citation Information

Xu Q, Niu L, Cao Z et al. **Genomic structural equation modelling of 12 cancers identifies a latent pan-cancer susceptibility factor and shared determinants of carcinogenesis** eBioMedicine, 2026; 129 (https://pubmed.ncbi.nlm.nih.gov/42361404/)

## Contact Information
Can Hou: houcan@wchscu.cn
