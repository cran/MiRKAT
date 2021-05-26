#' Microbiome Regression-based Analysis Test for a continuous outcome variable
#' 
#' Inner function for MiRKAT; computes MiRKAT for continuous outcomes. Called by MiRKAT if out_type="C"
#' 
#' This function is called by the exported function "MiRKAT" when the argument of MiRKAT, out_type, is set equal to "C". 
#' 
#' Each argument of MiRKAT_continuous is given the value of the corresponding argument given by the user to MiRKAT.
#' 
#' Function not exported 
#' 
#' @param y A numeric vector of the continuous outcome variable
#' @param X A numeric matrix or data frame containing additional covariates (default = NULL). If NULL, an intercept only model is used.
#' @param Ks A list of n by n kernel matrices (or a single n by n kernel matrix), where n is the sample size. It can be constructed 
#' from microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.
#' @param method A method to compute the kernel specific p-value (Default= "davies"). "davies" represents an exact method that computes
#'  the p-value by inverting the characteristic function of the mixture chisq. We adopt an exact variance component tests because most 
#'  of the studies concerning microbiome compositions have modest sample size. "moment" represents an approximation method that matches
#'  the first two moments. "permutation" represents a permutation approach for p-value calculation.
#' @param omnibus A string equal to either "Cauchy" or "permutation" (or nonambiguous abbreviations thereof), specifying whether 
#'  to use the Cauchy combination test or residual permutation to generate the omnibus p-value. 
#' @param nperm the number of permutations if method = "permutation" or when multiple kernels are considered. If method = "davies" or 
#' "moment", nperm is ignored. Defaults to 999.
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#'
#' @return If only one candidate kernel matrix is considered, returns a list containing the p-value for the candidate kernel matrix. 
#' If more than one candidate kernel matrix is considered, returns a list of two elements: the individual p-values for
#' each candidate kernel matrix, and the omnibus p-value
#'     \item{p_values}{p-value for each candidate kernel matrix}
#'     \item{omnibus_p}{omnibus p-value considering multiple candidate kernel matrices}
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'
#'
#'@author Ni Zhao
#'
#'@references 
#'Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J., Ringel, Y., Li, H. and Wu, M.C. (2015)).
#' Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807
#' 
#'Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and 
#'Microbiome Association Studies. 40: 5-19. doi: 10.1002/gepi.21934
#' 
#'Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal
#' Statistical Society. Series C , 29, 323-333.
#' 
#'Satterthwaite, F. (1946). An approximate distribution of estimates of variance components. Biom. Bull. 2, 110-114.
#' 
#'Lee S, Emond MJ, Bamshad MJ, Barnes KC, Rieder MJ, Nickerson DA; NHLBI GO Exome Sequencing Project-ESP Lung Project Team, 
#'Christiani DC, Wurfel MM, Lin X. (2012) Optimal unified approach for rare variant association testing with application to small
#' sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.
#' 
#'Zhou, J. J. and Zhou, H.(2015) Powerful Exact Variance Component Tests for the Small Sample Next Generation Sequencing Studies 
#'(eVCTest), in submission.


MiRKAT_continuous = function(y, X = NULL, Ks, method, omnibus, 
                             nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  
  om <- substring(tolower(omnibus), 1, 1)
  n = length(y) 
  
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
  
  KRVs = R2 = NULL 
  if (returnKRV | returnR2) {
    reskrv = scale(resid(lm(y ~ 1)))
    L = reskrv %*% t(reskrv)
    if (returnKRV) { KRVs <- unlist(lapply(Ks, FUN = function(k) calcKRVstat(k, L))) }
    if (returnR2) { R2 <- unlist(lapply(Ks, FUN = function(k) calcRsquared(k, L))) }
  }
  
  ## Prep RHS; take care of aliased variables and pivoting
  qX1 <- qr(X1)
  X1 <- X1[, qX1$pivot, drop=FALSE]
  X1 <- X1[, 1:qX1$rank, drop=FALSE]
  
  mod <- lm(y ~ X1-1)
  s2  = summary(mod)$s**2
  D0  = diag(n) 
  res = resid(mod)
  P0  = D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)
  px  = ncol(X1)
  
  ## Individual test statistics 
  Qs <- c() 
  for (i in 1:length(Ks)) {
    Qs <- c(Qs, getQ(K = Ks[[i]], res, s2))
  }
  
  ## Permuted test stats 
  if (method == "permutation" | (length(Ks) > 1 & om == "p")) {
    q_sim =  sapply(1:nperm, function(i) {
      ind <- sample(n)
      sapply(1:length(Ks),function(j){
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
      })
    })  
    q_sim = t(q_sim)
  }
  
  ## Individual p-values 
  if (method == "davies") {    
    lambda0 = lapply(Ks, getLambda_davies, P0)
    ps = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      ps[i] = getindivP_davies(Qs[[i]], lambda0[[i]], n, px)    
    }
    out_pvs = ps 
  } else if (method == "moment"){ 
    parm = sapply(Ks, getParamSat, P0, px)  
    p_sat = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      p_sat[i] = getSat(Qs[[i]], keppa_tlt = parm[,i]$keppa_tlt, niu_tlt = parm[,i]$niu_tlt)    
    }
    out_pvs = p_sat 
  } else if (method == substr("permutation", 1, nchar(method))) {
    # P = (# of statistics permutation >= observed test statistics + 1) / (#permutation + 1), 
    # thus the smallest p value will be  1/ (# permutation + 1).
    if (length(Ks) > 1) {
      Q_all = rbind(unlist(Qs), q_sim)
      p_all = 1 - (apply(Q_all, 2, rank) - 1)/(nperm + 1) 
      p_perm = p_all[1,]
      out_pvs = p_perm 
    } else {
      Q_all <- c(unlist(Qs), q_sim)
      p_all <- 1 - (rank(Q_all) - 1)/(nperm + 1)
      p_perm = p_all[1] 
      out_pvs <- p_perm
    }
    
  }
  
  ## name indiv p-values 
  names(out_pvs) = names(Ks) 
  
  if (length(out_pvs) == 1) {
    if (!returnKRV & !returnR2) {
      return(list(p_values = out_pvs))
    } else if (!returnKRV & returnR2) {
      return(list(p_values = out_pvs, R2 = R2))
    } else if (returnKRV & !returnR2) {
      return(list(p_values = out_pvs, KRV = KRVs))
    } else {
      return(list(p_values = out_pvs, KRV = KRVs, R2 = R2))    
    }
  }
  
  ## calculate omnibus p-value 
  if (om == "p") {
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1) 
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  
  } else if (om == "c") {
    cauchy.t <- sum(tan((0.5 - out_pvs)*pi))/length(out_pvs)
    p_final <- 1 - pcauchy(cauchy.t)
  } else {
    stop("I don't know that omnibus option. Please choose 'permutation' or 'Cauchy'.")
  }
  
  
  ## return all 
  if (is.null(KRVs) & is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final))
  } else if (is.null(KRVs) & !is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final, R2 = R2))
  } else if (!is.null(KRVs) & is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final, KRV = KRVs))
  } else {
    return(list(p_values = out_pvs , omnibus_p = p_final, KRV = KRVs, R2 = R2))    
  }
  
}
