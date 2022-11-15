#' Microbiome Regression-based Kernel Association Test
#'
#' Test for association between microbiome composition and a continuous or dichotomous outcome by incorporating 
#' phylogenetic or nonphylogenetic distance between different microbiomes.
#' 
#' y and X (if not NULL) should all be numeric matrices or vectors with the same number of rows.
#'
#' Ks should be a list of n by n matrices or a single matrix. If you have distance metric(s) from metagenomic data, each kernel can be
#' constructed through function D2K. Each kernel can also be constructed through other mathematical approaches.
#'
#' Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis
#'
#' Parameter "method" only concerns with how kernel specific p-values are generated. When Ks is a list of multiple kernels, omnibus 
#' p-value is computed through permutation from each individual p-values, which are calculated through method of choice.
#'
#' @param y A numeric vector of the a continuous or dichotomous outcome variable.
#' @param X A numeric matrix or data frame, containing additional covariates that you want to adjust for. If NULL, a intercept
#'  only model is used. Defaults to NULL.
#' @param Ks A list of n by n kernel matrices or a single n by n kernel matrix, where n is the sample size. It can be constructed 
#' from microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.
#' @param out_type An indicator of the outcome type ("C" for continuous, "D" for dichotomous). 
#' @param method Method used to compute the kernel specific p-value. "davies" represents an exact method that computes the p-value by
#'  inverting the characteristic function of the mixture chisq. We adopt an exact variance component tests because most of the studies
#'  concerning microbiome compositions have modest sample size. "moment" represents an approximation method that matches the first
#'  two moments. "permutation" represents a permutation approach for p-value calculation. Defaults to "davies".
#' @param omnibus A string equal to either "cauchy" or "permutation" (or nonambiguous abbreviations thereof), specifying whether 
#'  to use the Cauchy combination test or residual permutation to generate the omnibus p-value. 
#' @param nperm The number of permutations if method = "permutation" or when multiple kernels are considered. If method = "davies" or 
#' "moment", nperm is ignored. Defaults to 999.
#' @param returnKRV A logical indicating whether to return the KRV statistic (a measure of effect size). Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return R-squared. Defaults to FALSE. 
#'
#' @return 
#' Returns a list containing the following elements: 
#'    \item{p_values}{P-value for each candidate kernel matrix}
#'    \item{omnibus_p}{Omnibus p-value considering multiple candidate kernel matrices}
#'    \item{KRV}{Kernel RV statistic (a measure of effect size). Only returned if returnKRV = TRUE.}
#'    \item{R2}{R-squared. Only returned if returnR2 = TRUE.}
#'
#'@author 
#'Ni Zhao
#'
#'@references 
#'Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J., Ringel, Y., Li, H. and Wu, M.C. (2015)).
#' Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807
#'
#'Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome
#'Association Studies. 40: 5-19. doi: 10.1002/gepi.21934
#'
#'Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal 
#'Statistical Society. Series C , 29, 323-333.
#'
#'Satterthwaite, F. (1946). An approximate distribution of estimates of variance components. Biom. Bull. 2, 110-114.
#'
#'Lee S, Emond MJ, Bamshad MJ, Barnes KC, Rieder MJ, Nickerson DA; NHLBI GO Exome Sequencing Project-ESP Lung Project Team, 
#'Christiani DC, Wurfel MM, Lin X. (2012) Optimal unified approach for rare variant association testing with application to small 
#'sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.
#'
#'Zhou, J. J. and Zhou, H.(2015) Powerful Exact Variance Component Tests for the Small Sample Next Generation Sequencing Studies 
#'(eVCTest), in submission.
#'
#'@import MASS
#'@importFrom CompQuadForm davies
#'@importFrom GUniFrac GUniFrac
#'@import stats
#'
#'
#'@examples
#'library(GUniFrac)
#'
#'data(throat.tree)
#'data(throat.otu.tab)
#'data(throat.meta)
#'
#'unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
#'if (requireNamespace("vegan")) {
#'  library(vegan)
#'  BC= as.matrix(vegdist(throat.otu.tab, method="bray"))
#'  Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], BC = BC)
#'} else {
#'  Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"])
#'}
#'Ks = lapply(Ds, FUN = function(d) D2K(d))
#'
#'covar = cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male"))
#'
#'# Continuous phenotype
#'n = nrow(throat.meta)
#'y = rnorm(n)
#'MiRKAT(y, X = covar, Ks = Ks, out_type="C", method = "davies")
#'
#'# Binary phenotype 
#'y = as.numeric(runif(n) < 0.5)
#'MiRKAT(y, X = covar, Ks = Ks, out_type="D")
#'
#'
#' @export 
MiRKAT = function(y, X = NULL, Ks, out_type = "C", 
                  method = "davies", omnibus = "permutation", 
                  nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  
  method <- match.arg(method, choices = c("davies", "permutation", "moment"))
  omnibus <- match.arg(tolower(omnibus), choices = c("permutation", "cauchy"))
  
  n = length(y)
  
  if (any(is.na(y))){
    ids = which(is.na(y))
    stop(paste("missing response for subject(s)", ids, ", please remove before proceeding \n")) 
  }

  if(is.null(X)==FALSE){
    if(NROW(X)!= length(y)) stop("Dimensions of X and y don't match.")
  }
  
  if (is.matrix(Ks)) {
    Ks = list(Ks)
  }
  
  if (is.list(Ks)) {  
    if((any(lapply(Ks, "nrow")!= n))|(any(lapply(Ks,  "ncol")!= n))){
      stop("Distance matrix need to be n x n, where n is the sample size \n ")    
    } 
    if (!is.list(Ks)) {
      stop("Distance needs to be a list of n x n matrices or a single n x n matrix \n")  
    }
  }
  
  
  if (!is.null(X)){
    if (any(is.na(X))){
      stop("NAs in  covariates X, please impute or remove subjects with missing covariate values") 
    }  
  }
  
  if (method == "moment" & n < 100 & out_type == "C"){
    
    warning("Continuous outcome: sample size < 100, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  if (method == "moment" & n < 200 & out_type == "D"){
    
    warning("Continuous outcome: sample size < 200, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  
  if (!(out_type %in%  c("C", "D"))){
    stop("Only continuous and binary outcomes are supported by this function. Please choose out_type = \"C\" or \"D\", or use an alternative function for other outcome types.")
  }
  if(out_type  == "C"){
    re = MiRKAT_continuous(y, X = X, Ks = Ks, method = method, omnibus = omnibus, nperm = nperm, returnKRV = returnKRV, returnR2 = returnR2)  
  }
  
  if(out_type  == "D"){
    re = MiRKAT_binary(y, X = X, Ks = Ks, method = method, omnibus = omnibus, nperm = nperm, returnKRV = returnKRV, returnR2 = returnR2)  
  }
  
  return(re)
}

    