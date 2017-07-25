# Code obtained from Xiang Zhan at June/29/2017
# MMiRKAT and KRV both work for multivariate phenotypes. 
# However, MMiRKAT is recommended when the dimensionality of the phenotype is small; and KRV works for higher
# dimensional and/or structured phenotypes. 

### Required library
# library(CompQuadForm)

#########################################
### Simple auxiliary functions:
#########################################

## general the sqrt of a psd matrix
matsqrt=function(x, tol=10^-6){
  n=nrow(x)
  temp=eigen(x)
  val=temp$values
  vec=temp$vectors
  for(i in 1:n){ 
    if(abs(val[i])<tol) {val[i]=0
    } else val[i]=sqrt(val[i])
  }
  Dval=diag(val,n,n)
  out=vec%*%Dval%*%t(vec)
  return(out)
}

#########################################
### Major function:
#########################################


MMiRKAT <- function(Y, X=NULL, K) {
  if (any(is.na(Y))){
    ids = which(is.na(Y))
    stop(paste("subjects", ids, "has missing response, please remove before proceed \n")) 
  }
  
  if(is.null(X)==FALSE){
    if(NROW(X)!= NROW(Y)) stop("Dimensions of X and Y don't match.")
  }
  
  if (!is.null(X)){
    if (any(is.na(X))){
      stop("NAs in  covariates X, please impute or remove subjects which has missing covariates values") 
    }  
  }
  
  if(class(K) != "matrix") {
    stop("Please convert your kernel object into a matrix.")
  }
  
  if((nrow(K)!= nrow(Y))|(ncol(K)!= nrow(Y))){
    stop("kernel matrix need to be n x n, where n is the sample size \n ")    
  } 
  
  Sigma=cov(Y)
  iSig=solve(Sigma) ## requires n>p
  Y=Y%*%matsqrt(iSig)
  n <- ncol(K)
  if (is.null(X)) {
    X <- matrix(rep(1, n), n, 1)
  } else {
    X <- model.matrix(~., data=data.frame(X))
  }
  I <- diag(n)
  P <- I - X %*% solve(t(X) %*% X) %*% t(X)
  Y0 <- P %*% Y
  t0 <-  sum(Y0 %*% t(Y0) * K) / sum(Y0 * Y0)
  lambda1 <- eigen(P %*% K %*% P - t0 * P, only.values=T)$values
  lambda1 <- lambda1[abs(lambda1) >= 1e-10]
  lambda2 <- (svd(Y0/sqrt(n-1), nu=0, nv=0)$d)^2
  lambda2 <- lambda2[lambda2 >= 1e-10]
  lambdas <- as.vector(outer(lambda2, lambda1, '*'))
  pv=davies(0, lambda=sort(lambdas, decreasing = T))$Qq
  return(pv)
}

