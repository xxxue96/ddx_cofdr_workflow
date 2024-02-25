# ddx_cofdr_workflow
Differentiation and colocalization analysis at SNP- and gene-levels
## General description
1. Data preprocessing
2. Differentiation analysis
3. Colocalization analysis
## Working environment
python=3.7.12

R=4.2.1

Download scripts to your local working directory.
```
unzip ddx_cofdr_workflow.zip && cd ddx_cofdr_workflow
# Pleases modify paths of required files/tools or other paramaters in params.r and export_params.bash before any analysis.
```
## 1. Data preprocessing
Make subdirectories
```
mkdir -p 01_data/{00_Standardise_GWAS,01_LDSC_Results/{sumstats,h2,rg},02_tmvtnorm_Results,03_Gene_Result,04_FUMA_Results/{SNP2GENE,Cell_Type}}
```
Starndardise GWAS in R and save it to 01_data/00_Standardise_GWAS. 
```
source("01_data/01_Standardise_GWAS.r")
gwas = fread(raw_gwas_path) %>% mutate(FRQ = as.numeric(FRQ)) %>% select(SNP, CHR, BP, EFFECT_ALLELE, NON_EFFECT_ALLELE, FRQ, BETA, SE, P, Z, N, INFO)
gwas_formatted = Standardise_GWAS(gwas_sumstats_path=gwas, save_path=paste0("01_data/00_Standardise_GWAS/", trait, ".txt.gz"), ref_genome="GRCH37", convert_ref_genome=NULL, snp_ids_are_rs_ids = TRUE, INFO_filter = 0.3, FRQ_filter = 0.01, return_path = TRUE)
```
Run LDSC, and save rg and intercept with traits as rownames and colnames to 01_data/LDSC_Results. 
```
bash 04_script/02_01_Run_LDSC.bash t1 t2 # multiple traits are accepted
Rscript 04_script/02_02_Extract_LDSC.r t1 t2
```
Run mle.tmvnorm to estimate correlation coefficient due to sample overlap, and save estimates with traits as rownames and colnames to 01_data/02_tmvtnorm_Results. (For traits with invalidated LDSC estimates due to limited sample size)
```
Rscript 04_script/02_03_Fit_tmvtnorm.r t1 t2
```
Run MAGMA, TWAS, spTWAS, and PWAS for each trait in R, and save results to 01_data/03_Gene_Result.
```
rm(list=ls())
source("02_ddx/05_01_Run_DDx_Gene_MAGMA.r")
source("02_ddx/05_02_Run_DDx_Gene_PrediXcan.r")
source("02_ddx/05_03_Run_DDx_Gene_PWAS.r")
source("params.r") # please modify paramaters manually
FDRthres = 1 # keep all MAGMA/TWAS/spTWAS/PWAS-detected genes for further annotation
lapply(c(t1, t2), function(pheno){
	ddx_gwas = paste0("01_data/00_Standardise_GWAS/", pheno, ".txt.gz")
	twas_dir_ddx = paste0("01_data/03_Gene_Result/", pheno)
	magma_dir_ddx = paste0(twas_dir_ddx, "/01_MAGMA")
	pwas_dir_ddx = paste0(twas_dir_ddx, "/04_pwas")
	sapply(list(twas_dir_ddx, magma_dir_ddx, pwas_dir_ddx), dir.create, showWarnings = FALSE)

	Run_DDx_Gene_MAGMA(decor_gwas_paths=ddx_gwas, magma_dir=magma_dir_ddx, save_dir=twas_dir_ddx)
	Run_DDx_Gene_PrediXcan(traits=pheno, decor_gwas_paths=ddx_gwas, metaxcan_dir=twas_dir_ddx, save_dir=twas_dir_ddx)
	Run_DDx_Gene_PWAS(traits=pheno, decor_gwas_paths=ddx_gwas, pwas_dir=pwas_dir_ddx, save_dir=twas_dir_ddx)

	# remove temporary files
	file.remove(paste0(twas_dir_ddx,"/01_harmonized_gwas/",pheno,".txt.gz"))
})
```
Analyze gwas via FUMA "SNP2GENE" and "Cell Type" module. Save corresponding zipped results files to 01_data/04_FUMA_Results/{SNP2GENE,Cell_Type} and rename results files with trait names

## 2. Differeniation analysis
Make subdirectories
```
mkdir -p 02_ddx/{01_DDx_Results,02_DDx_SigSNP,03_DDx_Clump,04_FUMA_Results/{SNP2GENE,Cell_Type,Image/{01_SNP_Manhattan,02_Gene_Manhattan,03_MAGMA_Tissue,04_SNP_Annot,05_Cell_Type}},05_CCGWAS_Results,06_CCGWAS_SigSNP,07_SNP_Compare,08_DDx_Gene_Result,09_CCGWAS_Gene_Result,10_DDx_CCGWAS_Gene_Compare,11_Excel_Summary/{Figures,Tables}}
```
Run SNP- and gene-level differentiation analysis in R
```
source("02_ddx/DDx_Workflow.r")
source("02_ddx/DDx_Workflow_params.r")
# snp-level
ddx_wrap_up(t1, t2)

# Analyze ddx-derived gwas via FUMA, and save results to 02_ddx/04_FUMA_Results/{SNP2GENE,Cell_Type,Image}

# gene-level
ddx_wrap_up_gene(t1, t2)
```
Format results to readable excel files in R
```
source("02_ddx/06_01_DDx_Format_Excel.r")
source("02_ddx/06_02_DDx_Format_Image.r")
source("02_ddx/06_03_DDx_Excel_Summary.r")
source("02_ddx/06_03_DDx_Excel_Summary_more_info.r")
source("02_ddx/06_04_DDx_Count_Res.r")
# format differentiation results to excel
DDx_Format_Excel(t1, t2, save_dir="02_ddx/11_Excel_Summary/Tables")

# summarize all images from FUMA SNP2GENE and Cell Type module (multiple t1 and t2 are accepted)
DDx_Format_Image(t1, t2, save_path="02_ddx/11_Excel_Summary/Figures/all_trait_comparison_image.xlsx")

# summarize all excel files in excel_dir (multiple excel_dir are accepted)
DDx_Excel_Summary(excel_dir="02_ddx/11_Excel_Summary/Tables", save_path="02_ddx/11_Excel_Summary/all_trait_comparison_summary.xlsx")
DDx_Excel_Summary_more_info(excel_dir="02_ddx/11_Excel_Summary/Tables", save_path="02_ddx/11_Excel_Summary/all_trait_comparison_summary_more_info.xlsx")

# summarize number of significant results 
DDx_Count_Res(excel_dir="02_ddx/11_Excel_Summary/Tables", save_path="02_ddx/11_Excel_Summary/all_trait_comparison_DDx_Count.txt")
```
**Extension**: Multi-trait differentiation analysis

Make subdirectories
```
mkdir -p 02_ddx/{12_mtCOJO_Results,13_mtCOJO_Excel_Summary/{Figures,Tables}}
```
Run SNP- and gene-level differentiation analysis across multiple traits in R to find SNPs and genes specific to trait t1
```
source("02_ddx/mtCOJO_Workflow.r")
source("02_ddx/mtCOJO_Workflow_params.r")

# @mtcojo_name: set a pseudo-name for the mtCOJO-derived trait sumstats. Default is NULL to directly concatenate t1 and t2phenos as a pseudo-name.
mtCOJO_wrap_up(t1, t2phenos=c(t2,t3), mtcojo_name=NULL) # more t2phenos are accepted
mtCOJO_wrap_up_gene(t1, t2phenos=c(t2,t3), mtcojo_name=NULL)
```
Format mtCOJO results to readable excel files in R
```
source("02_ddx/07_01_mtcojo_Format_Excel.r")
source("02_ddx/07_02_mtcojo_Format_Image.r")
source("02_ddx/07_03_mtcojo_Excel_Summary.r")
source("02_ddx/07_03_mtcojo_Excel_Summary_more_info.r")
source("02_ddx/07_04_mtcojo_Count_Res.r")

t2phenos = c(t2, t3) # more t2phenos are accepted
traits = c(t1, t2phenos)

# format multi-trait differentiation results to excel
mtcojo_Format_Excel(t1, t2phenos, mtcojo_name=NULL, save_dir="13_mtCOJO_Excel_Summary/Tables")

# summarize all images from FUMA SNP2GENE and Cell Type module (directly use DDx_Format_Image())
DDx_Format_Image(t1, t2=t2phenos, save_path=paste0("02_ddx/13_mtCOJO_Excel_Summary/Figures/", t1, "_mtCOJO_image.xlsx"), ddx_name=NULL)

# summarize all excel files in excel_dir (multiple excel_dir are accepted)
mtCOJO_Excel_Summary(excel_dir="02_ddx/13_mtCOJO_Excel_Summary/Tables/", save_path=paste0("02_ddx/13_mtCOJO_Excel_Summary/", t1, "_mtCOJO_summary.xlsx"))
mtcojo_Excel_Summary_more_info(excel_dir="02_ddx/13_mtCOJO_Excel_Summary/Tables/", save_path=paste0("02_ddx/13_mtCOJO_Excel_Summary/", t1, "_mtCOJO_summary_more_info.xlsx")

# summarize number of significant results
mtcojo_Count_Res(excel_dir="02_ddx/13_mtCOJO_Excel_Summary/Tables", save_path=paste0("02_ddx/13_mtCOJO_Excel_Summary/Tables/",t1,"_mtCOJO_Count.txt"))
```
## 3. Colocalization analysis
Make subdirectories
```
mkdir -p 03_cofdr/{01_Cofdr_Result,02_Cofdr_SigSNP,03_Cofdr_Clump,04_Cofdr_Gene_Result,05_01_Cofdr_FUMA_Input,05_02_Cofdr_FUMA_Result/{SNP2GENE,GENE2FUNC},06_gwas_pw_Result,07_gwas_pw_SigSNP,08_HyPrColoc_Result,09_sigSNP_Comparison,10_Excel_Summary/{Figures,Tables},11_Multi_Cofdr_Summary/Tables}
```
Run SNP- and gene-level colocalization analysis
```
source("03_cofdr/Cofdr_Workflow.r")
source("03_cofdr/Cofdr_Workflow_params.r")
# snp-level
cofdr_wrap_up(t1, t2)

# Analyze cofdr-derived gwas in 03_cofdr/05_01_Cofdr_FUMA_Input via FUMA, and save results to 03_cofdr/05_02_Cofdr_FUMA_Result/{SNP2GENE,GENE2FUNC}

# gene-level
cofdr_wrap_up_gene(t1, t2)
```
Format results to readable excel files
```
source("03_cofdr/06_01_Cofdr_Format_Excel.r")
source("03_cofdr/06_03_Cofdr_Excel_Summary.r")
source("03_cofdr/06_03_Cofdr_Excel_Summary_more_info.r")
source("03_cofdr/06_04_Cofdr_Count_Res.r")
# format colocalization results to excel
Cofdr_Format_Excel(t1, t2, save_dir="03_cofdr/10_Excel_Summary/Tables")

# summarize all images from FUMA GENE2FUNC module (multiple t1 and t2 are accepted)
Cofdr_Format_Image(t1, t2, save_path = "03_cofdr/10_Excel_Summary/Figures/all_trait_comparison_image.xlsx")

# summarize all excel files in excel_dir (multiple excel_dir are accepted)
Cofdr_Excel_Summary(excel_dir="03_cofdr/10_Excel_Summary/Tables", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_summary.xlsx")
Cofdr_Excel_Summary_more_info(excel_dir="03_cofdr/10_Excel_Summary/Tables", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_summary_more_info.xlsx")

# summarize number of significant results
Cofdr_Count_Res(excel_dir="03_cofdr/10_Excel_Summary/Tables", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_Cofdr_Count.txt")
```






