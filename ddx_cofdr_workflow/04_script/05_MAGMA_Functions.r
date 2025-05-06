# aggregate snp-level effect to gene-level effect through positional mapping
library(MungeSumstats)
library(MAGMA.Celltyping)
library(data.table)
setDTthreads(40)
source("params.r")

# step1: ensure magma installation path
# old_path = Sys.getenv("PATH")
# Sys.setenv(PATH = paste(old_path, "/mnt/data/xue/Tool/MAGMA", sep=":"))
	
# step2: create magma storage file path	
get_magma_paths <- function(path=NULL, output_path=NULL, upstream_kb=0, downstream_kb=0, id_type="Entrez")
{
	gwasFileName = gsub("\\.txt.gz$", "", basename(path))
	pathMagmaFiles = sprintf("%s/%s.%sUP.%sDOWN/%s", output_path, gwasFileName, upstream_kb, downstream_kb, id_type)
	prefix = sprintf("%s.%sUP.%sDOWN", gwasFileName, upstream_kb, downstream_kb)
	filePathPrefix = sprintf("%s/%s.%sUP.%sDOWN", pathMagmaFiles, gwasFileName, upstream_kb, downstream_kb)
	dir.create(pathMagmaFiles, showWarnings = FALSE, recursive = TRUE)

	magmaPaths = list(files = pathMagmaFiles, prefix = prefix, filePathPrefix = filePathPrefix, gwasFileName = gwasFileName, gwasFilePath = path)
	
	return(magmaPaths)
}

# step3: gene-based analysis
# aggregate snp effects into a gene to detect disease-associated-genes
map_snps_to_genes <- function(path, genome_build=NULL, upstream_kb=0, downstream_kb=0, N=NULL, genome_ref_path=NULL, population="eur", genes_only=FALSE, storage_dir=tools::R_user_dir("MAGMA.Celltyping", which="cache"), force_new=FALSE, verbose=TRUE,
							  magmaPaths, genomeLocFile=NULL, genesAnnotPrefix=NULL)
{	
	# download 1000 Genomes reference panel if it's null
	if(is.null(genome_ref_path)){
		genome_ref_path = MAGMA.Celltyping::get_genome_ref(genome_ref_path = genome_ref_path, storage_dir = storage_dir, population = population, verbose = verbose)
	} else{
		print(paste("Use pre-downloaded 1kg reference panel:", basename(genome_ref_path)))
	}	
	# download NCBI gene loc file if it's null	
	if(is.null(genomeLocFile)){
		if(is.null(genome_build)){
			genome_build = MungeSumstats::get_genome_builds(sumstats_list = path, names_from_paths = TRUE)
		}
		if(toupper(genome_build) %in% c("GRCH36")){
			genomeLocFile = MAGMA.Celltyping:::get_genomeLocFile(build = "GRCH36", storage_dir = storage_dir)
		} else if(toupper(genome_build) %in% c("GRCH37","HG37","HG19")){
			genomeLocFile = MAGMA.Celltyping:::get_genomeLocFile(build = "GRCH37", storage_dir = storage_dir)
		} else if(toupper(genome_build) %in% c("GRCH38","HG38")){
			genomeLocFile = MAGMA.Celltyping:::get_genomeLocFile(build = "GRCH38", storage_dir = storage_dir)
		} else{
			stop("Genome build must be: 'GRCH36', `GRCH37', or 'GRCH38'")
		}
	} else{
		print(paste("Use pre-downloaded gene loc file:", basename(genomeLocFile)))
	}	
	# remove a trailing slash to avoid errors on windows
    outPath = gsub("\\/$", "", magmaPaths$filePathPrefix)
    genes_annot = sprintf("%s.genes.annot", genesAnnotPrefix)
    genes_out = sprintf("%s.genes.out", outPath)
	
	if((file.exists(genes_annot) & file.exists(genes_out)) & (force_new==FALSE)){
		message("Precomputed file detected: ", genes_out)
        return(genes_out)
	}
	# MAGMA requires files to be decompressed
    path = MAGMA.Celltyping:::decompress(path_formatted = path, remove = FALSE, overwrite = TRUE)
	# replcae empty column with 0.9999 in case of inconsistent ncol in different rows(may cause error in selecting SNP,P to create genes.out) https://stackoverflow.com/questions/38818613/how-to-execute-linux-commands-from-r-via-bash-under-the-windows-subsystem-for-li
	format_cmd = paste("sed -i -e 's/^\t/0.9999\t/' -e ':a' -e 's/\t\t/\t0.9999\t/g' -e 'ta' -e 's/\t$/\t0.9999/'", path)
	system(format_cmd)
	# check whether there is an N column in the sumstats file (if it wasn't provided as an argument)
	if(is.null(N)){
		first_line = readLines(path, n=1)
		column_headers = strsplit(first_line, "\t")[[1]]
		if("N" %in% column_headers){
			n_arg = "ncol=N"
		} else{
			nval = as.numeric(readline(paste("There is no N column within the sumstats file. What is the N value for this GWAS?")))
			if(is.na(nval)){stop(paste(nval, "provided but value of N for the GWAS must be numeric"))}
			if(nval < 1000){stop(paste("Value of N provided is less than 1,000. This seems unlikely."))}
			if(nval > 1e+08){stop(paste("Value of N provided is over than 00,000,000. This seems unlikely."))}
			n_arg = sprintf("N=%s", nval)
		}
	} else{
		n_arg = sprintf("N=%s", N)
	}	
    # create genes.annot
	if(!file.exists(genes_annot))
	{
		print(paste("Create genes.annot"))
		magma_cmd = sprintf(paste("magma", "--annotate window=%s,%s", "--snp-loc '%s'", "--gene-loc '%s'", "--out '%s'"),
							upstream_kb, downstream_kb,  paste0(genome_ref_path, ".bim"), genomeLocFile, genesAnnotPrefix)
		system(magma_cmd)
	}
	# create genes.out
	if(!file.exists(genes_out)){
		print(paste("Create genes.out"))
		magma_cmd = sprintf(paste("magma", "--bfile '%s'", "--pval '%s' use=SNP,P %s", if(isTRUE(genes_only)) "--genes-only" else NULL, "--gene-annot '%s.genes.annot'", "--out '%s'"),
							genome_ref_path, path, n_arg, genesAnnotPrefix, magmaPaths$filePathPrefix)	
		system(magma_cmd)
	}	
	# return path to genes.out file
	return(genes_out)
}

# step4: gene-set analysis
# aggregate gene effects for relationship between gene sets(binary phenotypes) and disease-associated-genes to detect disease-associated-gene-sets
map_genes_to_sets <- function(magmaPaths, upstream_kb=0, downstream_kb=0, set_name=NULL, set_path=NULL, force_new=FALSE)
{
	pathMagmaSet = sprintf("%s/%s", magmaPaths$files, set_name)
	dir.create(pathMagmaSet, showWarnings = FALSE, recursive = TRUE)
	
	outPath = sprintf("%s/%s", pathMagmaSet, magmaPaths$prefix)
    genes_raw = sprintf("%s.genes.raw", magmaPaths$filePathPrefix)
	sets_out = sprintf("%s.gsa.out", outPath)
	
	if((file.exists(genes_raw) & file.exists(sets_out)) & (force_new==FALSE)){
		message("Precomputed file detected: ", sets_out)
        return(sets_out)
	}
	# create gas.out for gene sets
	if(!file.exists(sets_out)){
		print(paste("Create gas.out for gene sets"))
		magma_cmd = paste("magma --gene-results", genes_raw, "--set-annot", set_path, "--out", outPath)
		system(magma_cmd)
	}
	# return path to gsa.out for gene sets
	return(sets_out)
}

# step5: gene property analysis (specificity analysis)
# tissue-specificity; aggregate gene effects for relationship between tissue specific gene expression profiles(continous phenotypes) and disease-associated-genes to detect disease-associated-tissue
# cell-type-specificity: scRNA data was downloaded from https://github.com/Kyoko-wtnb/FUMA_scRNA_data
map_genes_to_property <- function(magmaPaths, upstream_kb=0, downstream_kb=0, covar_name=NULL, covar_path=NULL, force_new=FALSE)
{	
	pathMagmaCovar = sprintf("%s/%s", magmaPaths$files, covar_name)
	dir.create(pathMagmaCovar, showWarnings = FALSE, recursive = TRUE)
	
	outPath = sprintf("%s/%s", pathMagmaCovar, magmaPaths$prefix)
    genes_raw = sprintf("%s.genes.raw", magmaPaths$filePathPrefix)
	covar_out = sprintf("%s.gsa.out", outPath)

	if((file.exists(genes_raw) & file.exists(covar_out)) & (force_new==FALSE)){
		message("Precomputed file detected: ", covar_out)
        return(covar_out)
	}
	# MAGMA requires files to be decompressed
    covar_path = MAGMA.Celltyping:::decompress(path_formatted = covar_path, remove = FALSE, overwrite = FALSE)
	# create gas.out for covariates (tissue expression)
	if(!file.exists(covar_out)){
		print(paste("Create gas.out for covariates"))
		magma_cmd = paste("magma --gene-results", genes_raw, "--gene-covar", covar_path, "--model condition-hide=Average direction=greater --out", outPath)
		system(magma_cmd)
	}
	# return path to gsa.out for covariates
	return(covar_out)
}

