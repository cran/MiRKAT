#' Kernel RV Coefficient Test; Inner Function
#' 
#' Function called when user calls function KRV. For each kernel matrix inputted into KRV, KRV runs inner.KRV on that 
#' kernel with the inputted kernel.y outcome matrix.
#' 
#'  y and X (if not NULL) should all be numerical matrices or vectors with the same number of rows.
#'
#'  Ks should be a list of n by n matrices or a single matrix. If you have distance metric from metagenomic data, each kernel can be
#'  constructed through function D2K. Each kernel can also be constructed through other mathematical approaches.
#'
#'  Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis
#'  
#' Parameter "method" only concerns how kernel specific p-values are generated. When Ks is a list of multiple kernels, omnibus
#' p-value is computed through permutation from each individual p-value, which are calculated through method of choice.
#' 
#' @param y  A numeric n by p matrix of p continuous phenotype variables and 
#' sample size n (default = NULL). If it is NULL, a 
#' phenotype kernel matrix must be entered for "kernel.y". Defaults to NULL.
#' @param X A numeric n by q matrix, containing q additional covariates 
#' (default = NULL). If NULL, an intercept only model is used. If the first
#' column of X is not uniformly 1, then an intercept column will be added. 
#' @param  adjust.type Possible values are "none" (default if X is null), 
#' "phenotype" to adjust only the y variable (only possible if y is a numeric 
#' phenotype matrix rather than a pre-computed kernel), or "both" to adjust 
#' both the X and Y kernels. 
#' @param kernel.otu A numeric OTU n by n kernel matrix or a list of matrices, 
#' where n is the sample size. It can be constructed from microbiome data, such 
#' as by transforming from a distance metric.
#' @param kernel.y Either a numerical n by n kernel matrix for phenotypes or a 
#' method to compute the kernel of phenotype. Methods are "Gaussian" or "linear". 
#' A Gaussian kernel (kernel.y="Gaussian") can capture the general relationship 
#' between microbiome and phenotypes; a linear kernel (kernel.y="linear") 
#' may be preferred if the underlying relationship is close to linear. 
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.
#'  
#' 
#' @return Returns a p-value for the candidate kernel matrix
#'     \item{pv}{p-value for the candidate kernel matrix}
#'     \item{KRV}{KRV statistic for the candidate kernel matrix. Only returned if returnKRV = TRUE.}  
#'     \item{R2}{R-squared for the candidate kernel matrix. Only returned if returnR2 = TRUE.}
#'     
#'@author 
#'Haotian Zheng, Xiang Zhan, Ni Zhao
#'
#'@references 
#' Zhan, X., Plantinga, A., Zhao, N., and Wu, M.C. A Fast Small-Sample Kernel Independence Test for Microbiome Community-Level 
#' Association Analysis. Biometrics. 2017 Mar 10. doi: 10.1111/biom.12684.
#' 
#'
#'
inner.KRV <- function(y = NULL, X = NULL, adjust.type, 
                      kernel.otu, kernel.y, returnKRV = FALSE, returnR2 = FALSE){
  
  K = kernel.otu 
  
  n=nrow(K)
  I.n=diag(1,n)
  I.1=rep(1,n)
  
  if (is.matrix(kernel.y)) {
    L = kernel.y
    if (!is.null(adjust.type)) {
      if (adjust.type == "b") {
        Px = X %*% chol2inv(chol(t(X) %*% X)) %*% t(X) # projection matrix of X
        K = (I.n - Px) %*% kernel.otu %*% (I.n - Px)
        L = (I.n - Px) %*% L %*% (I.n - Px)
      }
    }  
  } else {
    if (is.null(adjust.type)) {
      err.Y = y
      if (kernel.y == "Gaussian") { L = kern_g(err.Y) } 
      if(kernel.y == "linear") { L = err.Y%*%t(err.Y) }
    } else if (adjust.type == "b") {
      if (kernel.y == "Gaussian") { L = kern_g(y) } 
      if (kernel.y == "linear") { L = y %*% t(y) }
      Px = X %*% chol2inv(chol(t(X) %*% X)) %*% t(X) # projection matrix of X
      K = (I.n - Px) %*% kernel.otu %*% (I.n - Px)
      L = (I.n - Px) %*% L %*% (I.n - Px)
    } else if (adjust.type == "p") {
      Px=X%*%solve(t(X)%*%X)%*%t(X)
      err.Y=(I.n-Px)%*%y
      if (kernel.y == "Gaussian") { L = kern_g(err.Y) } 
      if(kernel.y == "linear") { L = err.Y%*%t(err.Y) }
    }
  }

  
  H=I.n-I.1%*%t(I.1)/n
  K=H%*%K%*%H
  L=H%*%L%*%H
  A=K/tr(K%*%K)  ## standard-version of K
  W=L/tr(L%*%L)
  
  Fstar=tr(A%*%W)
  mean.krv=tr(A)*tr(W)/(n-1)	## mean of KRV 
  
  T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
  Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
  temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
  temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
  temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
  temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
  temp2=temp21*temp22/temp23
  variance.krv=temp1+temp2		## variance of KRV
  
  T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
  T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
  t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
  t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
  t3=24*(n^2-n-4)*(U*Bs+B*Us)
  t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
  t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
  t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
  t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
  t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
  t8=24*(t81+t82)
  t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
  t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
  t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
  t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
  t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
  t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
  t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
  t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
  t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
  t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
  t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
  t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
  t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
  t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
  t20=-(n-2)*(t201+t202+t203)
  temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
  temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
  mom3=temp31/temp32
  skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
  
  m1=mean.krv
  m2=variance.krv
  m3=skewness.krv
  shape=4/m3^2
  scale=sqrt(m2)*m3/2
  location=m1-2*sqrt(m2)/m3
  PIIIpars=list(shape,location,scale)
  pv=1-ppearsonIII(Fstar, params=PIIIpars) 
  
  if (returnR2) {R2 = calcRsquared(K, L)} else {R2 = NULL}
  return(list(pv = pv, KRV = Fstar, R2 = R2))
}




## A toy function to calculate a Gaussian kernel matrix
kern_g = function(zz){
  n=nrow(zz)
  D=matrix(nrow=n, ncol=n)  ## the pairwise distance matrix
  for(i in 1:n){
    for(j in 1:n){
      D[i,j]=sum((zz[i,]-zz[j,])^2)
    }}	
  temp=c(D)
  D1=temp[temp>0]
  scl=median(D1)   ## use the median distance as the bandwidth
  K=matrix(nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      K[i,j]=exp(-sum((zz[i,]-zz[j,])^2)/scl)
    }}		
  return(K)
}




