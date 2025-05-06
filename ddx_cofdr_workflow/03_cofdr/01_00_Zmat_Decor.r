# remove impact of correlation due to sample overlap on Zscore matrix https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4859-7#Sec13
# @zmat: a n*d matrix of original z-scores with d equals to the number of SNPs common to n gwas studies; 
# @cor.mat: a n*n correlation matrix 

library(powerplus)

decor <- function(zmat, cor.mat)
{
	zdecor = Matpow(cor.mat,-0.5) %*% zmat		
	dimnames(zdecor) <- dimnames(zmat)
	
	return(zdecor)	
}
