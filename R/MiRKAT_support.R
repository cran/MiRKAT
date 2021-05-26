

# KRV test statistic 
calcKRVstat <- function(K, L) {
  n = nrow(K) 
  I.n=diag(1,n)
  I.1=rep(1,n)
  H=I.n-I.1%*%t(I.1)/n
  K=H%*%K%*%H
  L=H%*%L%*%H
  A=K/tr(K%*%K)  ## standard-version of K
  W=L/tr(L%*%L)
  Fstar=tr(A%*%W)
  return(Fstar)
}

# R-squared 
calcRsquared <- function(K, L) {
  r1 <- cor(as.numeric(K), as.numeric(L)) 
  return(r1^2)
}

# trace of a matrix 
tr=function(x){return(sum(diag(x))) }

# permutation matrix 
getPermuteMatrix <- function(perm, N, strata = NULL) {
  if (length(perm) == 1) {
    perm <- how(nperm = perm)
  }
  if (!missing(strata) && !is.null(strata)) {
    if (inherits(perm, "how") && is.null(getBlocks(perm))) 
      setBlocks(perm) <- strata
  }
  if (inherits(perm, "how")) 
    perm <- shuffleSet(N, control = perm)
  else {
    if (!is.integer(perm) && !all(perm == round(perm))) 
      stop("permutation matrix must be strictly integers: use round()")
  }
  if (is.null(attr(perm, "control"))) 
    attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
                                            nperm = nrow(perm)), class = "how")
  perm
}


# For binary outcome 
getHm = function(Q,muQ, varQ, df){
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p = 1 - pchisq(Q_corrected ,df = df)
  p = ifelse(p<0, 0, p)
  return(p)
}
Get_Var_Elements =function(m4,u1,u2){
  temp1 = u1^2 * u2^2
  a1    = sum(m4 * temp1)
  a2    = sum(u1^2) * sum(u2^2) - sum(temp1)
  a3    = sum(u1*u2)^2 - sum(temp1)
  return(a1+a2+2*a3)
}


