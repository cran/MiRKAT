#' Inner Function for CSKAT, Correlated Sequence Kernel Association Test
#'
#' Small-sample SKAT for correlated (continuous) data ('c' stands for 'correlated'). Computes the adjusted score statistic and p-value. 
#' 
#' @param formula.H0 A two-sided linear formula object under the null, indicating the variables to adjust. Only handles the mixed effects.
#'  Needed only if model = "gaussian" and method = "davies". 
#' @param data  an optional data frame containing the variables named in formula.
#' @param K    the kernel matrix, which quantifies the similarities between samples
#' @return
#' \describe{
#'   \item{p.value}{association p-value}
#'   \item{Q.adj}{adjusted score statistic}
#' }
#' @references
#' Zhan X, et al. (2018) A small-sample kernel association test for correlated data with application to microbiome association studies.  
#' Genet Epidemiol., submitted.
#' 
#' 


inner.CSKAT <- function (formula.H0, data = NULL, K) {
  
  H0.lmer <- formula.H0
  temp <- as.character(H0.lmer)
  temp[3] <- gsub('\\+\\s*\\(.*?\\)', '',  temp[3], perl=TRUE)
  
  H0.lm <- as.formula(paste(temp[2], temp[1], temp[3]))
  
  lmer.obj <- lmer(H0.lmer, data)
  
  var.d <- crossprod(as.matrix(getME(lmer.obj,"Lambdat")))
  Zt <- as.matrix(getME(lmer.obj, "Zt"))
  vr <- sigma(lmer.obj)^2
  
  var.b <- (t(Zt) %*% var.d %*% Zt)
  sI <- Diagonal(nrow(var.b))
  V2 <- vr * (var.b + sI)
  
  Vi <- sqrt.inv(V2)$Vi
  
  X1 <- model.matrix(H0.lm, data)
  lhs <- H0.lm[[2]]
  Y <- eval(lhs, data)
  
  X1 <- Vi %*% X1
  Y <- Vi %*% Y
  K <- Vi %*% K %*% Vi
  
  mod <- lm(Y ~ X1 - 1)
  res <- resid(mod)
  s2 <- sum(res^2) 
  
  D0  <- diag(length(Y)) 
  P0  <- D0 - X1 %*% solve(t(X1)%*%X1) %*% t(X1)
  PKP <- P0 %*% K %*% P0
  
  q <- as.numeric(res %*% K %*% res / s2)
  ee <- eigen(PKP - q * P0, symmetric = T)  
  lambda0 = ee$values[abs(ee$values) >= 1e-10]
  #	p1 <- davies(0, lambda=lambda0)$Qq
  p.value <- KAT.pval(0, lambda=sort(lambda0, decreasing=T))
  
  return(list(p.value = p.value, Q.adj = q, V = V2, resid = res))
}

