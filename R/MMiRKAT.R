#' Multivariate Microbiome Regression-based Kernel Association Test
#' 
#' Test for association between overall microbiome composition and multiple continuous outcomes. 
#'
#' Missing data is not permitted. Please remove all individuals with missing Y, X, K prior to analysis
#'
#' The method of generating kernel specific p-values is "davies", which represents an exact method that computes the p-value by
#' inverting the characteristic function of the mixture chisq.
#' 
#' @param Y A numerical n by p matrix of p continuous outcome variables, n being sample size.
#' @param X A numerical n by q matrix or data frame, containing q additional covariates that you want to adjust for (Default = NULL).
#' If it is NULL, an intercept only model is fit.
#' @param Ks A list of numerical n by n kernel matrices, or a single n by n kernel matrix, where n is the sample size. Kernels can be 
#' constructed from distance matrices (such as Bray-Curtis or UniFrac distances) using the function D2K, or through other mathematical approaches.  
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#' 
#' 
#' @return 
#' Returns a list of the MMiRKAT p-values for each inputted kernel matrix, labeled with the names of the kernels, if given.
#'     \item{p_values}{list of the p-values for each individual kernel matrix inputted}
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'
#' @author 
#' Nehemiah Wilson, Haotian Zheng, Xiang Zhan, Ni Zhao
#'
#' @references 
#' Zheng, H., Zhan, X., Tong, X., Zhao, N., Maity,A., Wu, M.C., and Chen,J. A small-sample multivariate kernel machine test for 
#' microbiome association studies. Genetic Epidemiology, 41(3), 210-220. DOI: 10.1002/gepi.22030 
#' 
#' 
#'@importFrom mixtools matsqrt
#'@examples  
#'library(GUniFrac)
#'if(requireNamespace("vegan")) { library(vegan) }
#'
#'data(throat.tree)
#'data(throat.otu.tab)
#'data(throat.meta)
#'
#'unifracs <- GUniFrac(throat.otu.tab, throat.tree, alpha=c(0, 0.5, 1))$unifracs
#'
#'if(requireNamespace("vegan")) {
#'  BC= as.matrix(vegdist(throat.otu.tab , method="bray"))
#'  Ds = list(w = unifracs[,,"d_1"], u = unifracs[,,"d_UW"], BC = BC) 
#'} else {
#'  Ds = list(w = unifracs[,,"d_1"], u = unifracs[,,"d_UW"]) 
#'}
#'Ks = lapply(Ds, FUN = function(d) D2K(d)) 
#'
#'n = nrow(throat.otu.tab)
#'Y = matrix(rnorm(n*3, 0, 1), n, 3)
#'
#'covar = cbind(as.numeric(throat.meta$Sex == "Male"), as.numeric(throat.meta$PackYears))
#'MMiRKAT(Y = Y, X = covar, Ks = Ks) 
#'
#'@export
#'
MMiRKAT <- function(Y, X = NULL, Ks, returnKRV = FALSE, returnR2 = FALSE){

  if (any(is.na(Y))){
    ids = which(is.na(Y))
    stop(paste("Subject(s)", ids, "has missing response, please remove before proceeding. \n")) 
  }
  
  if(is.null(X)==FALSE){
    if(NROW(X)!= NROW(Y)) stop("Dimensions of X and Y don't match.")
  }
  
  if (!is.null(X)){
    if (any(is.na(X))){
      stop("NAs in covariates X, please impute or remove subjects with missing covariates values.") 
    }  
  }
  
  if(!is.list(Ks)){
    Ks <- list(Ks)
  }
  
  if(!length(Ks)==1){
    if(is.null(names(Ks))){
      message("Your p-values are not labeled with their corresponding kernel matrix. In order to have them labeled,
              make your list of kernel matrices for the input of the form 'list(name1=K1, name2=K2'...) in order for the output
              p-values to be labeled with 'name1,' 'name2,' etc. \n ")
    }
  }
  
  lapply(Ks, FUN = function(K) {
    if (!is.matrix(K))  {
      stop("Please ensure that your kernel(s) are matrices.")
    }
    
    if ((nrow(K)!= nrow(Y)) | (ncol(K)!= nrow(Y))){
      stop("Dimension mismatch. Kernel matrix needs to be n x n, where n is the sample size. \n ")    
    } 
  })
  
  if (returnKRV | returnR2) {
    L1 = scale(Y) 
    L = L1 %*% t(L1) 
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
  
  pvals <- c()
  for(i in 1:length(Ks)){
    pvalue <- inner.MMiRKAT(Y = Y, X = X, K = Ks[[i]])
    pvals[i] <- pvalue
  }
  
  # Naming outputted p-values with the names of their corresponding kernel matrices
  names_Ks <- names(Ks)
  named_pvals <- setNames(pvals, names_Ks)
  
  if (!returnKRV & !returnR2) {
    return(list(p_values = named_pvals))
  } else if (!returnKRV & returnR2) {
    return(list(p_values = named_pvals, R2 = R2)) 
  } else if (returnKRV & !returnR2) {
    return(list(p_values = named_pvals, KRV = KRVs))
  } else {
    return(list(p_values = named_pvals, KRV = KRVs, R2 = R2))
  }
}
