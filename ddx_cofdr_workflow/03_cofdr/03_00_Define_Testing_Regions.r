# defining testing regions https://github.com/jrs95/hyprcoloc/issues/5

# method1: LD-defined regions 
## chunking up a large genomic region into multiple distinct LD-blocks https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/
# ldetect_data = paste0("/mnt/data/xue/Tool/gwas-pw-0.21/ldetect_data/fourier_ls-all_", pop, ".bed")


# method2: COVID-associated loci-defined regions when performing colocalization COVID and the other trait
## https://www.sciencedirect.com/science/article/pii/S0197458022001117


# method3: drug-traget/eqtl defined regions when perform colocalization between qtl and the other trait