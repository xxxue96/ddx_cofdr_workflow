# @N: number of analyzed traits
# @zmat: merged zmat with SNP as rows and z-score for each phenotype as columna
cofdr <- function(N,zmat, iter=500, tol=1e-8, empNULL = FALSE, nulltype=1, bre = 120, df = 7){

  #***************************
  # any number of traits
  #***************************

  ##install packages (to be done once)
  #install.packages("locfdr")
  #install.packages("devtools")
  # library(devtools, quietly=TRUE)
  #install_github('jgscott/FDRreg', subdir="R_pkg/")

  library(FDRreg, quietly=TRUE)
  library(locfdr, quietly=TRUE)
  
  totalsnps= nrow(zmat)
  empnull.mean = rep(0,N)
  empnull.sd = rep(1,N)
  if (empNULL==TRUE){
    ##empirical null estimation
    for (i in 1:N) {
      locobj.za = locfdr(zmat[,i], nulltype=nulltype, bre = bre, df = df) 
      empnull.mean[i] = locobj.za$fp0[3,1]
      empnull.sd[i] = locobj.za$fp0[3,2]
    }
  }
  
  
  #no. of H0 and H1 combinations = 2^N
  k = 2^N
  
  totalSNPs=totalsnps
  f0z.mat = matrix(nrow=totalSNPs, ncol=N)  ##each column is f0_z for each trait
  f1z.mat = matrix(nrow=totalSNPs, ncol=N)  ##each column is f1_z for each trait
  pi0 = vector()  # vector of estimated pi0 for each trait
  
  for (i in 1:N) {
    probj.z = prfdr(zmat[,i], mu0=empnull.mean[i], sig0=empnull.sd[i]) ## theoretical null, can fit empirical null by locfdr 
    f0z.mat[,i]= probj.z$f0_z
    f1z.mat[,i]= probj.z$f1_z
    pi0[i]=probj.z$pi0
    #probj.zb = prfdr(zb, mu0=0, sig0=1) ## theoretical null, can fit empirical null by locfdr 
  }
  
  list_N = rep(  list(c(0,1)),N   )
  grid = expand.grid(list_N)  ##should have 2^N = k rows
  grid= grid+1
  #     Var1 Var2 Var3
  #1    1    1    1
  #2    2    1    1
  #3    1    2    1
  #4    2    2    1
  #5    1    1    2
  #6    2    1    2
  #7    1    2    2
  #8    2    2    2
  
  
  list.f01z = list(f0z.mat=f0z.mat, f1z.mat=f1z.mat)
  pi1 = 1-pi0
  list.pi0  = list(pi0=pi0,pi1=pi1)
  
  #fj = f(joint distribution)
  fjoint=matrix(nrow=totalSNPs, ncol=k)
  for (i in 1:k){
    
    prod1=rep(1,totalSNPs)
    for  (j in 1:N) {
      grid.entry = grid[i,j]
      prod1 = prod1* list.f01z[[grid.entry]][,j]    ##grid.entry=1 -->ie choose f0z, if entry=2, choose f1z
    }
    fjoint[,i] = prod1
  }
  #fjoint: each column is the joint distribution of eg f000,f001,f010 etc
  ##f00, f01, etc. will remain unchanged in EM iterations
  
  
  pi.joint = matrix(nrow=iter+1,ncol=k)
  #pi.joint = matrix()
  
  for (i in 1:k){
    prod1=1
    for  (j in 1:N) {
      grid.entry = grid[i,j]
      prod1 = prod1*list.pi0[[grid.entry]][j]     ##grid.entry=1 -->ie choose f0z, if entry=2, choose f1z
    }
    pi.joint[1,i] = prod1
  }
  ##pi.joint = each row is each EM iteration; each column is pi00, pi01, p10, pi11 etc (posterior probability of each combination of H0 and h1)
  
  nosnps = totalSNPs
  diff.pi.all0 = 1e5
  #********************************
  # EM loop
  #********************************
  
  
  t1 = proc.time()
  
  m=1
  while ((m<=iter)&&(diff.pi.all0>tol) ) {
    
    #initilize
    Emat = matrix(nrow=totalSNPs, ncol=k)
    
    ##E-step 
    f.zazb = 0
    for (i in 1:k){
      f.zazb = f.zazb + pi.joint[m,i]*fjoint[,i]
      #f.zazb = pi00[m]*f00 + pi01[m]*f01 + pi10[m]*f10 + pi11[m]*f11
    }
    
    #for (p in 1:nosnps) {  
    for (i in 1:k){
      #Emat[p,i]= pi.joint[m,i]*fjoint[,i][p]/f.zazb[p]
      Emat[,i]= pi.joint[m,i]*fjoint[,i]/f.zazb
    } 
    #}
    
    ## M step 
    pi.vector = apply(Emat,2,mean)
    m <- m+1
    for (i in 1:k) {
      pi.joint[m,i]=pi.vector[i]
    }
    #pi01[m] = pi.vector[2]
    #pi10[m] = pi.vector[3]
    #pi11[m] = pi.vector[4]
    
    #diff.pi.all0 = abs(pi.joint[m,1] - pi.joint[m-1,1])  
    #use sum of pi differences as tolerance measure (ie sum up all differences of pi000,pi100,pi001 etc)
    diff.pi.all0 = sum( abs(pi.joint[m,] - pi.joint[m-1,])  ) 
    #if (m%%50==0) {cat(m,"iterations ...",'\n')}
    #print(m)
  }
  
  #proc.time()-t1
  
  #print(pi.vector)
  
  ##check convergence
  #print(diff.pi.all0)
  
  
  
  #***************************************************
  # likehood ratio test 
  # H0: pi11 = pi1* times pi*1
  # H0 
  #**************************************************
  #pi under H0 
  
  ind=NULL
  pi.star = NULL
  ## find those
  for (i in 1:N){
    ##look at each column of grid, find out elements starting with a "2" in the grid 
    
    pi.star[i] = sum(pi.vector[which(grid[,i]==2)])
  }
  
  piAll_1.H0 <- prod(pi.star)  
  
  
  #log likelihood under H1
  pi.vector.mult <- function(row) {sum(pi.vector*row)}
  snp.likelihood = apply(fjoint, 1, pi.vector.mult )  
  likeli.H1 = sum( log(snp.likelihood))  
  
  
  #log likelihood under H0
  pi.vector.reserved = pi.vector[-c(1,length(pi.vector))]  ##reserve all elements except the 1st and the last ie pi000 and pi111
  pi.vector.H0 = c( 1-piAll_1.H0-sum(pi.vector.reserved) , pi.vector.reserved,  piAll_1.H0 )
  
  pi.vector.mult.H0 <- function(row) {sum(pi.vector.H0*row)}
  snp.likelihood = apply(fjoint, 1, pi.vector.mult.H0)  
  likeli.H0 = sum( log(snp.likelihood)  )  
  
  ##likelihood ratio statistic D
  D = 2*(likeli.H1-likeli.H0) 
  LRT.p = pchisq(D, df=1, ncp = 0, lower.tail = FALSE)
  
  #column names of posterior prob. matrix 
  colname =NULL
  
  for (r in 1:nrow(grid)) {
    colname[r] <- paste("comb",grid[r,]-1,collapse="")
  }
  colnames(Emat) = colname
  return( list(pi.vector=pi.vector, post.prob=Emat, diff = diff.pi.all0, D=D, LRT.p = LRT.p  )   )
}
