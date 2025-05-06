#********************************************************************
# Infer the summary statistics for differential diagnoses(ddx) based
# on summary statistics of involved traits
# @param:ldsc_intercept(LD score regression intercept)
# @param: gwas1(a dataframe stored the summary statistics of trait 1)
# @param: t1.case(sample size of trait 1)
# @param: gwas2(a dataframe stored the summary statistics of trait 2)
# @param: t2.case(sample size of trait 2)
# return: ddx.gwas(summary statistics for ddx)
# Note: the dataframe for summary statistics should include beta(OR),
# se, eaf,snpid, chr, bpos, a1, a2 columns
#********************************************************************
sddx <- function(ldsc_intercept, gwas1, t1.case, gwas2, t2.case){
  #obtain beta and se from 2 sets of gwas summary statistics
  if ('BETA' %in% colnames(gwas1)){beta1=gwas1$BETA} else if 
     ('OR' %in% colnames(gwas1)){beta1=log(gwas1$OR)} else {
      cat('please check input gwas1 file: no beta and OR column',sep='\n');break
  }
  if ('SE' %in% colnames(gwas1)){se1=gwas1$SE} else {
  cat('please check input gwas1 file: no se column',sep='\n');break
  }
      
  if ('BETA' %in% colnames(gwas2)){beta2=gwas2$BETA} else if 
     ('OR' %in% colnames(gwas2)){beta2=log(gwas2$OR)} else {
      cat('please check input gwas2 file: no beta and OR column',sep='\n');break
  }
  if ('SE' %in% colnames(gwas2)){se2=gwas2$SE} else {
  cat('please check input gwas2 file: no se column',sep='\n');break
  }
  
  var.beta1=se1^2
  var.beta2=se2^2
  
  #obtain estimated cov.beta12 from LSDC program
  cor.beta12.est=ldsc_intercept
  cov.beta12=cor.beta12.est*se1*se2
    
  #calculate differential beta, se, zval and pval between 2 sets of gwas summary statistics
  ddx_beta=beta1-beta2
  ddx_se=sqrt(var.beta1+var.beta2-2*cov.beta12)
  ddx_zval=ddx_beta/ddx_se
  ddx_pval=pnorm(-abs(ddx_zval))*2
  
  #calculate sample size and effect allele frequency from 2 sets of original gwas datasets
  ddx.n=t1.case+t2.case
  if ('FRQ' %in% colnames(gwas1)&'FRQ' %in% colnames(gwas2)) {
  ddx.eaf=(gwas1$FRQ*gwas1$N+gwas2$FRQ*gwas2$N)/(gwas1$N+gwas2$N)  
  
  attach(gwas1)
  ddx.gwas <- data.frame(SNP=SNP,CHR=CHR,BP=BP,A1=A1,A2=A2,FRQ=ddx.eaf,BETA=ddx_beta,SE=ddx_se,P=ddx_pval,Z=ddx_beta/ddx_se,N=ddx.n) 
  detach(gwas1)} else {
  
  cat('there is no eaf column in gwas1 or gwas2',sep='\n')
  attach(gwas1)
  ddx.gwas <- data.frame(SNP=SNP,CHR=CHR,BP=BP,A1=A1,A2=A2,BETA=ddx_beta,SE=ddx_se,P=ddx_pval,Z=ddx_beta/ddx_se,N=ddx.n) 
  detach(gwas1)
  }
  return(ddx.gwas)
}