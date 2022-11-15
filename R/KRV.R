#' Kernel RV Coefficient Test (KRV)
#' 
#' Kernel RV coefficient test to evaluate the overall association between 
#' microbiome composition and high-dimensional or structured phenotype or genotype.
#' 
#' kernels.otu should be a list of numerical n by n kernel matrices, or a single 
#' n by n kernel matrix, where n is sample size.
#'
#' When kernel.y is a method ("Gaussian" or "linear") to compute the kernel of 
#' phenotype, y should be a numerical phenotype matrix, and X (if not NULL) 
#' should be a numeric matrix of covariates. Both y and X should have n rows.
#'
#' When kernel.y is a kernel matrix for the phenotype, there is no need to provide 
#' X and y, and they will be ignored if provided. In this case, kernel.y and 
#' kernel.otu should both be numeric matrices with the same number of rows and columns.
#'
#' Missing data is not permitted. Please remove all individuals with missing 
#' kernel.otu, y (if not NULL), X (if not NULL), and kernel.y (if a matrix is 
#' entered) prior to analysis.
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
#' @param kernels.otu A numeric OTU n by n kernel matrix or a list of matrices, 
#' where n is the sample size. It can be constructed from microbiome data, such 
#' as by transforming from a distance metric.
#' @param kernel.y Either a numerical n by n kernel matrix for phenotypes or a 
#' method to compute the kernel of phenotype. Methods are "Gaussian" or "linear". 
#' A Gaussian kernel (kernel.y="Gaussian") can capture the general relationship 
#' between microbiome and phenotypes; a linear kernel (kernel.y="linear") 
#' may be preferred if the underlying relationship is close to linear. 
#' @param omnibus A string equal to either "Cauchy" or "kernel_om" (or unambiguous 
#' abbreviations thereof), specifying whether to use the Cauchy combination test 
#' or an omnibus kernel to generate the omnibus p-value. 
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
#' Liu, Hongjiao, Ling, W., Hua, X., Moon, J.Y., Williams-Nguyen, J., Zhan, X., Plantinga, A.M., Zhao, N., 
#' Zhang, A., Durazo-Arzivu, R.A., Knight, R., Qi, Q., Burk, R.D., Kaplan, R.C., and Wu, M.C. 
#' Kernel-based genetic association analysis for microbiome phenotypes identifies host genetic 
#' drivers of beta-diversity. 2021+ 
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
#' ## Simulate covariate data 
#' set.seed(123)
#' n = nrow(throat.otu.tab)
#' Sex <- throat.meta$Sex
#' Smoker <- throat.meta$SmokingStatus
#' anti <- throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage
#' Male = (Sex == "Male")**2
#' Smoker = (Smoker == "Smoker") **2
#' Anti =  (anti != "None")^2
#' cova = cbind(1, Male, Smoker, Anti)
#'
#' ## Simulate microbiome data
#' otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff
#' unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs
#' # Distance matrices
#' D.weighted = unifracs[,,"d_1"]
#' D.unweighted = unifracs[,,"d_UW"]
#' # Kernel matrices 
#' K.weighted = D2K(D.weighted)
#' K.unweighted = D2K(D.unweighted)
#' 
#' if (requireNamespace("vegan")) {
#'  library(vegan)
#'  D.BC = as.matrix(vegdist(otu.tab.rff, method="bray"))
#'  K.BC = D2K(D.BC)
#'} 
#'
#' # Simulate phenotype data 
#' rho = 0.2
#' Va = matrix(rep(rho, (2*n)^2), 2*n, 2*n)+diag(1-rho, 2*n)
#' Phe = mvrnorm(n, rep(0, 2*n), Va)  
#' K.y = Phe %*% t(Phe)  # phenotype kernel 
#'
#' # Simulate genotype data
#' G = matrix(rbinom(n*10, 2, 0.1), n, 10) 
#' K.g = G %*% t(G)  # genotype kernel
#'
#'
#' ## Unadjusted analysis (microbiome and phenotype)
#' KRV(y = Phe, kernels.otu = K.weighted, kernel.y = "Gaussian")  # numeric y 
#' KRV(kernels.otu = K.weighted, kernel.y = K.y)  # kernel y 
#' 
#' ## Adjusted analysis (phenotype only)
#' KRV(kernels.otu = K.weighted, y = Phe, kernel.y = "linear", X = cova, adjust.type = "phenotype")
#' 
#' if (requireNamespace("vegan")) {
#'   ## Adjusted analysis (adjust both kernels; microbiome and phenotype) 
#'   KRV(kernels.otu = K.BC, kernel.y = K.y, X = cova, adjust.type='both')
#' 
#'   ## Adjusted analysis (adjust both kernels; microbiome and genotype)
#'   KRV(kernels.otu = K.BC, kernel.y = K.g, X = cova, adjust.type='both')
#' }
#'
#'
#'@export 
KRV <- function(y = NULL, X = NULL, adjust.type = NULL, kernels.otu, kernel.y, 
                omnibus = "kernel_om", returnKRV = FALSE, returnR2 = FALSE){
  
  ## Check input 
  
  om <- substring(tolower(omnibus), 1, 1)
  
  if (!is.null(adjust.type)) {
    adjust.type <- substring(tolower(adjust.type), 1, 1)
    
    if (adjust.type == "p") {
      if (!kernel.y %in% c("linear", "Gaussian")) {
        stop("Phenotype-only adjustment may only be used with numeric y, not pre-computed phenotype kernel. Try option 'both' instead or double check that you have entered either 'linear' or 'Gaussian' for kernel.y.")
      } 
    } 
    
    if (!adjust.type %in% c("p", "b")) {
      stop("I don't know that covariate adjustment choice. Please choose 'phenotype' or 'both', or leave this option NULL.")
    }
    
  } 
  
  if (!is.list(kernels.otu)) {
    kernels.otu <- list(kernels.otu)
  }
  
  if (!all(unlist(lapply(kernels.otu, FUN = function(k) is.matrix(k))))) {
    stop("Please ensure kernels.otu is either a single n x n kernel matrix or a list of n x n kernel matrices.")
  }
  
  if (!is.null(X) & !all(X[,1] == 1)) {
    X1 <- cbind(1, X) 
    X <- X1 
  }
  
  if (is.null(X) & !is.null(adjust.type)) {
    warning("X is NULL, so no covariate adjustment will be done.")
    adjust.type = NULL 
  }
  
  if (is.matrix(kernel.y)){
    n = nrow(kernels.otu[[1]])
    if (!is.null(y)){
      warning("When a phenotype kernel is provided, argument \"y\" will be ignored.\n")
    }
    if (ncol(kernels.otu[[1]]) != n | nrow(kernel.y) != n | ncol(kernel.y) != n){
      stop("Kernel matrices need to be n x n, where n is the sample size.\n ")
    }
  } else if (!is.matrix(kernel.y)){
    if (!(kernel.y %in%  c("Gaussian", "linear"))){
      stop("Please choose kernel.y = \"Gaussian\" or \"linear\", or enter a kernel matrix for \"kernel.y\".\n")
    }
    if(is.null(y)){
      stop("Please enter a phenotype matrix for argument \"y\" or enter a kernel matrix for argument \"kernel.y\".\n")
    }
    n = NROW(y)
    if (!all(unlist(lapply(kernels.otu, FUN = function(x) nrow(x) == n & ncol(x) == n)))) {
      stop("Kernel matrix/matrices must be n x n, where n is the sample size. \n ")
    }
    if (any(is.na(y))){
      ids = which(is.na(y))
      stop(paste("Missing response for subject(s)", ids, "- please remove before proceeding. \n"))
    }
    if (!is.null(X)){
      if (any(is.na(X))){
        stop("NAs in covariates X, please impute or remove subjects with missing covariates values.\n") 
      }  
      if(NROW(X)!= NROW(y)) stop("Dimensions of X and y don't match.\n")
    }
  }
  
  
  ## Actual test 
  
  pvals <- c()
  if (returnKRV) { KRVs <- c() } else { KRVs = NULL }
  if (returnR2) { R2 <- c() } else { R2 = NULL }
  
  for(i in 1:length(kernels.otu)){
    res <- inner.KRV(y = y, X = X, adjust.type = adjust.type, kernel.otu=kernels.otu[[i]], 
                     kernel.y = kernel.y, returnKRV = TRUE, returnR2 = TRUE)
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
    if (om == "k") {
      K.om <- matrix(0, nrow = nrow(kernels.otu[[1]]), ncol = ncol(kernels.otu[[1]]))
      for(i in 1:length(kernels.otu)){
        K.om = K.om + kernels.otu[[i]]/tr(kernels.otu[[i]])
      }
      omnibus_p <- as.numeric(inner.KRV(y = y, X = X, adjust.type = adjust.type, 
                                        kernel.otu = K.om, kernel.y = kernel.y)$pv)
    } else if (om == "c") {
      cauchy.t <- sum(tan((0.5 - pvals)*pi))/length(pvals)
      omnibus_p <- 1 - pcauchy(cauchy.t)
    } else {
      stop("I don't know that omnibus option. Please choose 'kernel_om' or 'Cauchy'.")
    }
    
    ## Return for multiple kernels
    if (is.null(KRVs) & is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = omnibus_p))
    } else if (is.null(KRVs) & !is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = omnibus_p, R2 = R2))
    } else if (!is.null(KRVs) & is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = omnibus_p, KRV = KRVs))
    } else {
      return(list(p_values = pvals, omnibus_p = omnibus_p, KRV = KRVs, R2 = R2))    
    }
  }

  
  ## Return for single kernels
  if (is.null(KRVs) & is.null(R2)) {
    return(list(p_values = pvals))
  } else if (is.null(KRVs) & !is.null(R2)) {
    return(list(p_values = pvals, R2 = R2))
  } else if (!is.null(KRVs) & is.null(R2)) {
    return(list(p_values = pvals, KRV = KRVs))
  } else {
    return(list(p_values = pvals, KRV = KRVs, R2 = R2))    
  }
  
}

