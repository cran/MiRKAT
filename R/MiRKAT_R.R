#'Robust MiRKAT (robust regression)
#'
#'A more robust version of MiRKAT utilizing a linear model by robust regression using an M estimator. 
#'
#' MiRKAT.R creates a kernel matrix using the linear model created with the function rlm, a robust regression function, then does 
#' the KRV analysis on Ks and the newly formed kernel matrix representing the outome traits. 
#'
#' y and X should all be numerical matrices or vectors with the same number of rows, and mustn't be NULL.
#'
#' Ks should be a list of n by n matrices or a single matrix. If you have distance metric from metagenomic data, each kernel can be 
#' constructed through function D2K. Each kernel may also be constructed through other mathematical approaches.
#'
#' Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis
#'
#'@param y A numeric vector of the a continuous or dichotomous outcome variable.
#'@param X A numerical matrix or data frame, containing additional covariates that you want to adjust for Mustn't be NULL
#'@param Ks list of n by n kernel matrices (or a single n by n kernel matrix), where n is the sample size. It can be constructed from 
#'microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#'@return
#'Returns p-values for each individual kernel matrix, an omnibus p-value if multiple kernels were provided, and measures of effect size KRV and R2. 
#'    \item{p_values}{labeled individual p-values for each kernel}
#'    \item{omnibus_p}{omnibus p_value, calculated as for the KRV test}
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'     
#'@author
#'Weijia Fu
#'
#'
#'
#'
#'@examples     
#'
#'
#'# Generate data
#'library(GUniFrac)
#'data(throat.tree)
#'data(throat.otu.tab)
#'data(throat.meta)
#'
#'unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
#'Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], 
#'          BC= as.matrix(vegdist(throat.otu.tab, method="bray"))) 
#'Ks = lapply(Ds, FUN = function(d) D2K(d))
#'
#'covar = cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male"))
#'
#'# Continuous phenotype
#'n = nrow(throat.meta)
#'y = rchisq(n, 2)
#'MiRKAT.R(y, X = covar, Ks = Ks)
#'
#'@export
#'
MiRKAT.R <- function(y, X, Ks, returnKRV = FALSE, returnR2 = FALSE){
  if (is.null(X)) {
    stop("Please provide a covariate matrix X. To fit an intercept-only model, use MiRKAT instead of MiRKAT-R, as no robust 
         intercept-only model is available.") 
  } 
  
  if (returnKRV | returnR2) {
    warning("Note that R2 and KRV are calculated using an intercept-only model, and therefore will be identical between MiRKAT-R and MiRKAT.")
    reskrv = scale(resid(lm(y ~ 1)))
    L = reskrv %*% t(reskrv)
    if (returnKRV) {
      KRVs <- unlist(lapply(Ks, FUN = function(k) calcKRVstat(k, L)))
    } else {
      KRVs = NULL 
    }
    if (returnR2) {
      R2 <- unlist(lapply(Ks, FUN = function(k) calcRsquared(k, L)))
    } else {
      R2 = NULL 
    }
  } else {
    KRVs = R2 = NULL 
  }
  

  mod <- rlm(y ~ X)
  res <- resid(mod)
  scalep<- mod$s
  k <- 1.345
  res0 <- res/scalep
  u <- ifelse(abs(res0)>k, k*res0/abs(res0), res0)
  U <- u %*% t(u)
  sig <- KRV(kernels.otu=Ks, kernel.y=U)
  
  return(list(p_values = sig$p_values, omnibus_p = sig$omnibus_p, KRV = KRVs, R2 = R2))
}

