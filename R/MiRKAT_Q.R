#'Robust MiRKAT (quantile regression)
#'
#'A more robust version of MiRKAT utlizing a linear model that uses quantile regression. 
#'
#'MiRKAT.Q creates a kernel matrix using the linear model created with the function rq, a quantile regression function, then does 
#'the KRV analysis on Ks and the newly formed kernel matrix representing the outome traits. 
#'
#'Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis
#'
#'@param y A numeric vector of the a continuous or dichotomous outcome variable.
#'@param X A numerical matrix or data frame, containing additional covariates that you want to adjust for. Mustn't be NULL.
#'@param Ks A list of n by n kernel matrices (or a single n by n kernel matrix), where n is the sample size. If you have distance metric from metagenomic data, each kernel can 
#' be constructed through function D2K. Each kernel can also be constructed through other mathematical approaches, such as linear or Gaussian kernels.
#'@param omnibus A string equal to either "Cauchy" or "kernel_om" (or nonambiguous abbreviations thereof), specifying whether 
#'  to use the Cauchy combination test or an omnibus kernel to generate the omnibus p-value. 
#'@param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#'@param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#' 
#'@return
#'Returns p-values for each individual kernel matrix, an omnibus p-value if multiple kernels were provided, and measures of effect size KRV and R2. 
#'    \item{p_values}{labeled individual p-values for each kernel}
#'    \item{omnibus_p}{omnibus p_value, calculated as for the KRV test}
#'    \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'    \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'
#'@author
#'Weija Fu
#'
#'@importFrom quantreg rq
#'
#'@examples     
#'library(quantreg)
#'library(vegan) 
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
#'covar = scale(cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male")))
#'
#'# Continuous phenotype
#'n = nrow(throat.meta)
#'y = rchisq(n, 2) + apply(covar, 1, sum) 
#'MiRKAT.Q(y, X = covar, Ks = Ks)
#'
#'@export
MiRKAT.Q <-function(y, X, Ks, omnibus = "kernel_om", returnKRV = FALSE, returnR2 = FALSE){
  
  om <- substring(tolower(omnibus), 1, 1)
  
  if (is.null(X)) {
    stop("Please provide a covariate matrix X to use robust MiRKAT. To fit an intercept-only model, use MiRKAT instead of MiRKAT-R, 
as no robust intercept-only model is available. \n") 
  }
  
  if (is.matrix(Ks)) {
    Ks = list(Ks)
  }
  
  if (returnKRV | returnR2) {
    warning("Note that R2 and KRV are calculated using an intercept-only model, and therefore will be identical between MiRKAT-Q and MiRKAT \n.")
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
  
  mod <- rq(y~X)
  res <- resid(mod)
  res <- ifelse(res<0, 0.5, -0.5)
  res <- as.matrix(res)
  U <- res %*% t(res)
  sig <- KRV(kernels.otu = Ks, kernel.y=U, omnibus = omnibus)
  
  
  if (length(Ks) == 1) {
    if (!returnKRV & !returnR2) {
      return(list(p_values = sig$p_values))
    } else if (!returnKRV & returnR2) {
      return(list(p_values = sig$p_values, R2 = R2))
    } else if (returnKRV & !returnR2) {
      return(list(p_values = sig$p_values, KRV = KRVs))
    } else {
      return(list(p_values = sig$p_values, KRV = KRVs, R2 = R2))
    }
  }
  if (!returnKRV & !returnR2) {
    return(list(p_values = sig$p_values, omnibus_p = sig$omnibus_p))
  } else if (!returnKRV & returnR2) {
    return(list(p_values = sig$p_values, omnibus_p = sig$omnibus_p, R2 = R2))
  } else if (returnKRV & !returnR2) {
    return(list(p_values = sig$p_values, omnibus_p = sig$omnibus_p, KRV = KRVs))
  } else {
    return(list(p_values = sig$p_values, omnibus_p = sig$omnibus_p, KRV = KRVs, R2 = R2))
  }
}


