# Consensus-Based Gene Mapping

This folder contains the public workflow scripts for the consensus-based gene mapping analysis of pan-cancer GWAS risk loci. The workflow integrates eight evidence sources: eQTL colocalization, sQTL colocalization, blood pQTL colocalization, protein-altering variant linkage, ClinVar, MGI, H-MAGMA, and PoPS.

## Scripts

- `consensus_mapping_pipeline.sh`: connects the analysis steps with editable path placeholders. QTL data extraction and H-MAGMA parameter generation are not included; the script starts from prepared analysis inputs. The H-MAGMA block is retained as a one-tissue template and can be repeated across annotation tissues.
- `coloc_eqtl.R`, `coloc_sqtl.R`, `coloc_pqtl.R`: run coloc on prepared locus-level GWAS-QTL harmonized files.
- `pav.py`, `clump_merge.py`, `pav_select.py`: identify protein-altering variants in LD with lead or independent significant variants.
- `clinvar.py`, `mgi.py`: query ClinVar and MGI phenotype resources for cancer-relevant evidence near pan-cancer loci.
- `HMAGMA_merge.py`: merge MAGMA gene-level outputs generated with H-MAGMA annotations.
- `pops_merge.py`: merge PoPS predictions and retain locus-level gene ranks.
- `consensus.py`: aggregate all evidence sources and write consensus gene-mapping outputs.
- `consensus_dict_pan.txt`: tissue, phenotype, and annotation groups used for the pan-cancer consensus score.

## Expected Inputs

The workflow expects `loci.parameter` from the verified pan-cancer loci, a harmonized pan-cancer GWAS file with `snp` and `p` columns, prepared locus-level QTL colocalization input files, QTL parameter files, a European LD reference panel, VEP annotations, ClinVar and MGI database files, H-MAGMA annotation files, and standard PoPS resources.

Prepared QTL parameter files should contain whitespace-delimited rows defining `prefix`, `tissue`, locus or locus list, GWAS sample size, QTL sample size, and trait type. The pQTL parameter file additionally includes the data-preparation script and QTL data path fields retained from the original workflow.
