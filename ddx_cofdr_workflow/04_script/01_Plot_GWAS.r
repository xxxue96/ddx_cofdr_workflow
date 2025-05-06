library(qqman)
library(data.table)
setDTthreads(40)

Plot_GWAS <_ function(gwas_sumstats_path)
{
	jpeg(file = "QQ_Manhattan_Plot.jpeg")
	par(mfrow = c(1, 2))

	gwas = fread(gwas_sumstats_path)
	qq(gwas$P, main = "QQ_Plot", col = "blue4", pch = 19, cex = 0.6, las = 1)
	manhattan(gwas, main = "Manhattan_Plot", col = c("blue4", "orange3"), cex = 0.6, cex.axis = 0.9)
	
	# Calculating Genomic Inflation Factor
	z = gwas$Z
	lambda = round(median(z^2) / 0.454, 3)
	print(paste("Genomic Inflation Factor:", lambda))
	
	dev.off()
}

