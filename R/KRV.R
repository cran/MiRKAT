#' Kernel RV Coefficient Test (KRV)
#' 
#' Kernel RV coefficient test to evaluate the overall association between microbiome composition and high-dimensional or
#' structured phenotype.
#' 
#' kernels.otu should be a list of numerical n by n kernel matrices, or a single n by n kernel matrix, where n is sample size.
#'
#' When kernel.y is a method ("Gaussian" or "linear") to compute the kernel of phenotype, y should be a numerical phenotype matrix, 
#' and X (if not NULL) should be a numeric matrix of covariates. Both y and X should have n rows.
#'
#' When kernel.y is a kernel matrix for the phenotype, there is no need to provide X and y, and they will be ignored if provided.
#' In this case, kernel.y and kernel.otu should both be numeric matrices with the same number of rows and columns.
#'
#' Missing data is not permitted. Please remove all individuals with missing kernel.otu, y (if not NULL), X (if not NULL), and 
#' kernel.y (if a matrix is entered) prior to analysis.
#'
#' @param y  A numeric n by p matrix of p continuous phenotype variables and sample size n (default = NULL). If it is NULL, a 
#' phenotype kernel matrix must be entered for "kernel.y". Defaults to NULL.
#' @param X A numeric n by q matrix, containing q additional covariates (default = NULL). If NULL, 
#'an intercept only model is used. No covariate adjustment is possible if a matrix is provided in kernel.y.  
#' @param kernels.otu A numeric OTU n by n kernel matrix or a list of matrices, where n is the sample size. It can be 
#' constructed from microbiome data, such as by transforming from a distance metric.
#' @param kernel.y Either a numerical n by n kernel matrix for phenotypes or a method to compute the kernel of phenotype. Methods are "Gaussian" or "linear". 
#' A Gaussian kernel (kernel.y="Gaussian") can capture the general relationship between microbiome and phenotypes; a linear kernel (kernel.y="linear") 
#' may be preferred if the underlying relationship is close to linear. 
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#' 
#' @return 
#' If only one candidate kernel matrix is considered, returns a list containing the p-value for the candidate kernel matrix. 
#' If more than one candidate kernel matrix is considered, returns a list of two elements: 
#'     \item{p_values}{P-value for each candidate kernel matrix}
#'     \item{omnibus_p}{Omnibus p-value} 
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'
#'@author
#'Nehemiah Wilson, Haotian Zheng, Xiang Zhan, Ni Zhao
#'
#'@references
#' Zheng, Haotian, Zhan, X., Plantinga, A., Zhao, N., and Wu, M.C. A Fast Small-Sample Kernel Independence Test for Microbiome 
#' Community-Level Association Analysis. Biometrics. 2017 Mar 10. doi: 10.1111/biom.12684.
#'
#'   
#'@importFrom PearsonDS ppearsonIII
#'
#'@examples         
#'library(GUniFrac)
#'library(MASS)
#'
#'data(throat.tree)
#'data(throat.otu.tab)
#'data(throat.meta)
#'
#'
#'set.seed(123)
#'n = nrow(throat.otu.tab)
#'Sex <- throat.meta$Sex
#'Smoker <- throat.meta$SmokingStatus
#'anti <- throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage
#'Male = (Sex == "Male")**2
#'Smoker =(Smoker == "Smoker") **2
#'anti =  (anti != "None")^2
#'cova = cbind(Male,  anti)
#'
#'otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff
#'unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs
#'
#'D.weighted = unifracs[,,"d_1"]
#'D.unweighted = unifracs[,,"d_UW"]
#'D.BC= as.matrix(vegdist(otu.tab.rff , method="bray"))
#'
#'
#'K.weighted = D2K(D.weighted)
#'K.unweighted = D2K(D.unweighted)
#'K.BC = D2K(D.BC)
#'
#'rho = 0.2
#'Va = matrix(rep(rho, (2*n)^2), 2*n, 2*n)+diag(1-rho, 2*n)
#'G = mvrnorm(n, rep(0, 2*n), Va)
#'
#'KRV(kernels.otu = K.weighted, kernel.y = G %*% t(G)) 
#'
#'
#'
#'@export 
KRV <- function(y = NULL, X = NULL, kernels.otu, kernel.y, returnKRV = FALSE, returnR2 = FALSE){
  
  if(!is.list(kernels.otu)){
    kernels.otu <- list(kernels.otu)
  }
  if(!length(kernels.otu)==1){
    if(is.null(names(kernels.otu))){
      message("The p-values will not be labeled according to their corresponding kernel matrices. For labeled p-values,
please form your list of kernels for input via 'list(name=K1, name=K2...)' to label the p-values with the corresponding names.")
    }
  }
  
  pvals <- c()
  if (returnKRV) { KRVs <- c() } else { KRVs = NULL }
  if (returnR2) { R2 <- c() } else { R2 = NULL }
  
  for(i in 1:length(kernels.otu)){
    res <- inner.KRV(y = y, X = X, kernel.otu=kernels.otu[[i]], kernel.y = kernel.y, returnKRV = TRUE, returnR2 = TRUE)
    pvals[i] <- res$pv
    if (returnKRV) { KRVs[i] <- res$KRV }
    if (returnR2) { R2[i] <- res$R2 } 
  }
  
  
  # Naming the p-values with the Kernel matrix names from Ks
  kernel.names <- names(kernels.otu)
  names(pvals) <- kernel.names
  if (returnKRV) { names(KRVs) <- kernel.names }
  if (returnR2) { names(R2) <- kernel.names }
  
  #Omnibus Test
  if (length(kernels.otu) > 1) {
    K.om <- matrix(0, nrow = nrow(kernels.otu[[1]]), ncol = ncol(kernels.otu[[1]]))
    for(i in 1:length(kernels.otu)){
      K.om = K.om + kernels.otu[[i]]/tr(kernels.otu[[i]])
    }
    
    omnibus_p <- as.numeric(inner.KRV(y = y, X = X, kernel.otu = K.om, kernel.y = kernel.y)$pv)
    return(list(p_values = pvals, omnibus_p = omnibus_p, KRV = KRVs, R2 = R2))
  }
  
  return(list(p_values = pvals, KRV = KRVs, R2 = R2))
  
}

