library(data.table)
library(stringr)
source("02_ddx/DDx_Workflow_params.r")

# step1: run metaxcan 
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

# step2: extract sig.genes
extract_gene_metaxcan <- function(traits, metaxcanPaths, save_dir)
{
	qtl_ann = fread(sqtl_ann)
	
	filter_genes <- function(genes_out_path)
	{
		tissue = str_split_fixed(basename(genes_out_path), "__PM__", 2)[, 2]
		if(tissue == ""){
			tissue = "smultixcan"
			dat = fread(genes_out_path)[, c("gene_name", "pvalue", "z_mean")]
			colnames(dat)[3] = "zscore"
		} else{
			tissue = gsub(".csv", "", tissue)
			dat = fread(genes_out_path)[, c("gene_name", "pvalue", "zscore")]
		}
		dat$tissue = tissue
		
		dat$FDR = p.adjust(dat$pvalue, method="fdr")
		dat = dat[dat$FDR<FDRthres,]
		dat = dat[order(dat$FDR, decreasing=FALSE),]			
		return(dat)
	}
	all_genes =  lapply(metaxcanPaths, function(x) lapply(x, filter_genes))
	
	## save eGenes
	sig_ddx_eqtl = do.call(rbind.data.frame, all_genes$predixcan_eqtl)
	sig_ddx_eqtl = rbind(sig_ddx_eqtl, all_genes$multixcan_eqtl[[1]])
	sig_ddx_eqtl = sig_ddx_eqtl[order(sig_ddx_eqtl$FDR, decreasing=FALSE),]
	colnames(sig_ddx_eqtl)[1] = "hgnc_symbol"
	path_sig_ddx_eqtl = paste0(save_dir, "/", traits, "_sig.eqtl.txt.gz")
	fwrite(sig_ddx_eqtl[, c("hgnc_symbol", "pvalue", "tissue", "FDR", "zscore")], file=path_sig_ddx_eqtl, sep="\t")
	
	# save and annotate sGenes
	sig_ddx_sqtl = do.call(rbind.data.frame, all_genes$predixcan_sqtl)
	sig_ddx_sqtl = rbind(sig_ddx_sqtl, all_genes$multixcan_sqtl[[1]])
	colnames(sig_ddx_sqtl)[1] = "sqtl_hg38"
	sig_ddx_sqtl = merge(sig_ddx_sqtl, qtl_ann, by="sqtl_hg38")
	sig_ddx_sqtl = sig_ddx_sqtl[order(sig_ddx_sqtl$FDR, decreasing=FALSE),]
	path_sig_ddx_sqtl = paste0(save_dir, "/", traits, "_sig.sqtl.txt.gz")
	fwrite(sig_ddx_sqtl[, c("sqtl_hg38", "pvalue", "tissue", "FDR", "hgnc_symbol", "sqtl_hg19", "ensemble_id", "zscore")], file=path_sig_ddx_sqtl, sep="\t")
	
	return(list(eqtl=path_sig_ddx_eqtl, sqtl=path_sig_ddx_sqtl))
}

Run_DDx_Gene_PrediXcan <- function(traits, decor_gwas_paths, metaxcan_dir, save_dir)
{
	metaxcanPaths = create_metaxcan_gene_results(traits, decor_gwas_paths, metaxcan_dir)
	ddx_gene = extract_gene_metaxcan(traits, metaxcanPaths, save_dir)
	
	return(ddx_gene)
}