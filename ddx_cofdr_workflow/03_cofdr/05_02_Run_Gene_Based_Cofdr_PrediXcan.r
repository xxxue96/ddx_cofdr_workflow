# colocalizaed eGenes/sGenes (i.e. map snps to genes based on eql/sqtl)
library(dplyr)
library(stringr)
library(data.table)
setDTthreads(40)
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_00_Zmat_Decor.r")
source("03_cofdr/01_02_Run_Cofdr.r")
source("03_cofdr/Cofdr_Workflow_params.r")

# step1: decorrelate gwas summary statistics
# @decor_path: directory to save decorrelated gwas, "A2" as effect allele
decor_snp4metaxcan <- function(traits, gwasFilePaths, merged_path=NULL, cor.mat=NULL, decor_path=tempdir())
{
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	if(!is.null(cor.mat)) {zmat = t(decor(t(zmat), cor.mat))}

	## save decorrelated gwas with SNP, P, N
	common_snp = common_snp[, c("SNP", "CHR", "BP", "A1", "A2", as.vector(outer(c("FRQ", "SE", "N"), traits, paste, sep="_"))), with=FALSE]
	common_snp = cbind(as.data.frame(common_snp), as.data.frame(zmat))		
	create_decor_gwas_paths <- function(pheno)
	{
		decor_gwas = common_snp[, c("SNP", "CHR", "BP", "A1", "A2", paste0(c("FRQ_", "SE_", "N_", "Z_"), pheno))]
		colnames(decor_gwas) = c("SNP", "CHR", "BP", "A1", "A2", "FRQ", "SE", "N", "Z")
		
		decor_gwas$BETA = decor_gwas$Z * decor_gwas$SE
		decor_gwas$P = 2*pnorm(-abs(decor_gwas$Z))
		
		outPath = paste0(decor_path, "/", pheno, ".txt.gz")
		fwrite(decor_gwas, file=outPath, sep="\t", quote=FALSE)
		return(outPath)
	}
	decor_gwas_paths = lapply(traits, create_decor_gwas_paths)
	
	return(decor_gwas_paths)
}


# step2: run metaxcan 
# @decor_gwas_paths: decorrelated gwas paths
# @metaxcan_dir: directory to save metaxcan results
create_metaxcan_gene_results <- function(traits, decor_gwas_paths, metaxcan_dir)
{
	run_predixcan <- function(i)
	{
		pheno=traits[i]; INPUT=dirname(decor_gwas_paths[[i]]); OUTPUT=metaxcan_dir; OUTPUT_Res=metaxcan_dir		
		cmd = paste("bash 04_script/06_02_Run_PrediXcan.bash", pheno, INPUT, OUTPUT, OUTPUT_Res)
		system(cmd)
		
		eqtl = list.files(path=paste0(metaxcan_dir, "/02_spredixcan/eqtl/", pheno), pattern=".csv", full.names=TRUE)
		sqtl = list.files(path=paste0(metaxcan_dir, "/02_spredixcan/sqtl/", pheno), pattern=".csv", full.names=TRUE)		
		return(list(eqtl, sqtl))
	}
	predixcan_res = lapply(1:length(traits), run_predixcan)
	predixcan_eqtl = do.call(c, lapply(predixcan_res, `[[`, 1))
	predixcan_sqtl = do.call(c, lapply(predixcan_res, `[[`, 2))
	
	multixcan_eqtl = list.files(path=paste0(metaxcan_dir, "/03_smultixcan/eqtl"), pattern=".txt$", full.names=TRUE)
	multixcan_sqtl = list.files(path=paste0(metaxcan_dir, "/03_smultixcan/sqtl"), pattern=".txt$", full.names=TRUE)	
	
	metaxcanPaths=list(predixcan_eqtl=predixcan_eqtl, predixcan_sqtl=predixcan_sqtl, multixcan_eqtl=multixcan_eqtl, multixcan_sqtl=multixcan_sqtl)
	return(metaxcanPaths)
}


# step3: run gene-based cofdr with zmat obtained from step2
# @metaxcanPaths: metaxcan result paths
# @save_dir: directory to save cofdr results	
run_metaxcan_cofdr <- function(traits, metaxcanPaths, save_dir)
{	
	pheno = paste(traits, collapse="_")
	qtl_ann = fread(sqtl_ann)
	
	## extract zscore for predixcan results
	extract_predixcan_z <- function(genes_out_path)
	{
		dat = fread(genes_out_path)[, c("gene_name", "zscore")]
		pheno = str_split_fixed(basename(genes_out_path), "__PM__", 2)[, 1]
		colnames(dat) = c("SNP", paste0("Z_",pheno))
			
		return(dat)
	}
	## cal. z-score from pvalue (metaxcan fle format https://predictdb.org/post/2022/03/08/metaxcan-output-file-formats/ ) for multixcan results
	convert_multixcan_p_to_z <- function(genes_out_path)
	{
		dat = fread(genes_out_path)[, c("gene", "pvalue")]
		pheno = gsub("_smultixcan.txt", "", basename(genes_out_path))		
		
		set.seed(nchar(pheno)*nchar(genes_out_path))
		dat$Z = abs(qnorm(dat$pvalue/2)) * sample(x=c(-1,1), size=length(dat$pvalue), replace=TRUE)
		
		dat = dat[, c("gene", "Z")]
		colnames(dat) = c("SNP", paste0("Z_",pheno))
			
		return(dat)
	}
	## run cofdr for metaxcan gene zscores, remove NA/Inf for zscore cols
	get_cofdr_res <- function(gwas_list, save_path)
	{
		merged_path = gwas_list %>% reduce(inner_join, by = "SNP")
		merged_path = merged_path[is.finite(rowSums(merged_path[, -"SNP"])),]
		gene_based_cofdr = Run_Cofdr(traits, gwasFilePaths=NULL, merged_path, cor.mat=NULL, save_path, keep_Zstat=TRUE)
			
		return(gene_based_cofdr$allsnp)
	}
	## run cofdr for predixcan gene zscores
	get_cofdr_res_predixcan <- function(m)
	{
		eqtl_path = grep(m, metaxcanPaths$predixcan_eqtl, value=TRUE)
		eqtl_res = get_cofdr_res(gwas_list=lapply(eqtl_path, extract_predixcan_z), save_path=paste0(save_dir, "/02_spredixcan/eqtl/", pheno, "_", m, ".txt.gz"))
		
		sqtl_path = grep(m, metaxcanPaths$predixcan_sqtl, value=TRUE)
		sqtl_res = get_cofdr_res(gwas_list=lapply(sqtl_path, extract_predixcan_z), save_path=paste0(save_dir, "/02_spredixcan/sqtl/", pheno, "_", m, ".txt.gz"))
		
		return(list(eqtl_res, sqtl_res))
	}
	## filter and annot sig.eqtl
	annot_eqtl <- function(cofdr_res, convert_gene_id=FALSE)
	{
		dat = fread(cofdr_res)[, c("SNP", "TT")]
		if(convert_gene_id == TRUE){colnames(dat)[1] = "ensemble_id"} else{colnames(dat)[1] = "hgnc_symbol"}
		
		dat = dat[dat$TT>=TTthres_eqtl, ]
		dat$tissue = qdapRegex::ex_between(basename(cofdr_res), paste0(pheno,"_"), ".txt.gz")[[1]]
		
		if(convert_gene_id == TRUE)
		{
			dat = merge(dat, qtl_ann[, c("hgnc_symbol", "ensemble_id")], by="ensemble_id")
			dat = dat[, -"ensemble_id"]
			dat = dat[!duplicated(dat), ]
		} 
			
		dat = dat[order(dat$TT, decreasing=TRUE),]
		return(dat[, c("hgnc_symbol", "TT", "tissue")])
	}
	## filter and annot sig.sqtl (hg38-based)
	annot_sqtl <- function(cofdr_res)
	{
		dat = fread(cofdr_res)[, c("SNP", "TT")]; colnames(dat)[1] = "sqtl_hg38"
		
		dat = dat[dat$TT>=TTthres_sqtl, ]
		dat$tissue = qdapRegex::ex_between(basename(cofdr_res), paste0(pheno,"_"), ".txt.gz")[[1]]
			
		dat = merge(dat, qtl_ann, by="sqtl_hg38")
		dat = dat[order(dat$TT, decreasing=TRUE),]
		return(dat)
	}	
	
	## multixcan across 49 tissues
	multixcanGWAS_eqtl = lapply(metaxcanPaths$multixcan_eqtl, convert_multixcan_p_to_z)
	multixcanGWAS_sqtl = lapply(metaxcanPaths$multixcan_sqtl, convert_multixcan_p_to_z)
	cofdr_eqtl = get_cofdr_res(gwas_list=multixcanGWAS_eqtl, save_path=paste0(save_dir, "/03_smultixcan/eqtl/", pheno, "_smultixcan.txt.gz"))
	cofdr_sqtl = get_cofdr_res(gwas_list=multixcanGWAS_sqtl, save_path=paste0(save_dir, "/03_smultixcan/sqtl/", pheno, "_smultixcan.txt.gz"))
	
	## predixcan for each tissue
	predixcan_cofdr_res = apply(model_src, 1, get_cofdr_res_predixcan)
	predixcan_cofdr_eqtl = do.call(c, lapply(predixcan_cofdr_res, `[[`, 1))
	predixcan_cofdr_sqtl = do.call(c, lapply(predixcan_cofdr_res, `[[`, 2))
		
	## filter cofdr_eqtl
	sig_cofdr_eqtl = do.call(rbind.data.frame, lapply(predixcan_cofdr_eqtl, annot_eqtl, convert_gene_id=FALSE))
	sig_cofdr_eqtl = rbind(annot_eqtl(cofdr_eqtl, convert_gene_id=TRUE), sig_cofdr_eqtl)	
	
	## filter and annotate cofdr_sqtl with annottaion file generated by leafcutter	
	sig_cofdr_sqtl = do.call(rbind.data.frame, lapply(c(cofdr_sqtl, predixcan_cofdr_sqtl), annot_sqtl))
	
	path_sig_cofdr_eqtl = paste0(save_dir, "/", pheno, "_sig.eqtl.txt.gz"); fwrite(sig_cofdr_eqtl, file=path_sig_cofdr_eqtl, sep="\t", quote=FALSE)
	path_sig_cofdr_sqtl = paste0(save_dir, "/", pheno, "_sig.sqtl.txt.gz"); fwrite(sig_cofdr_sqtl, file=path_sig_cofdr_sqtl, sep="\t", quote=FALSE)
	
	return(list(eqtl=path_sig_cofdr_eqtl, sqtl=path_sig_cofdr_sqtl))
}

# wrap-up functions
Run_Gene_Based_Cofdr_MetaXcan <- function(traits, gwasFilePaths, merged_path, cor.mat=NULL, save_dir)
{	
	decor_gwas_paths = decor_snp4metaxcan(traits, gwasFilePaths, merged_path, cor.mat, decor_path=tempdir())
	metaxcanPaths = create_metaxcan_gene_results(traits, decor_gwas_paths, save_dir)
	gene_based_cofdr = run_metaxcan_cofdr(traits, metaxcanPaths, save_dir)
	
	return(gene_based_cofdr)
}