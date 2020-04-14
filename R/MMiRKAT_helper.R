# Function to calculate one MMiRKAT p-value 
inner.MMiRKAT <- function(Y, X, K) {
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
  pv=CompQuadForm::davies(0, lambda=sort(lambdas, decreasing = TRUE))$Qq
  return(pv)
}

