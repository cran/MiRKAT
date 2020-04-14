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
#' @param nperm the number of permutations if method = "permutation" or when multiple kernels are considered. If method = "davies" or 
#' "moment", nperm is ignored. Defaults to 999.
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#'
#' 
#' @return If only one candidate kernel matrix is considered, returns a list containing the p-value for the candidate kernel matrix. 
#' If more than one candidate kernel matrix is considered, returns a list of two elements: the individual p-values for
#' each candidate kernel matrix, and the omnibus p-value
#'     \item{indivP}{p-value for each candidate kernel matrix}
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


MiRKAT_continuous = function(y, X = NULL, Ks, method = "davies", nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  n = length(y) 
  
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
  
  if (returnKRV | returnR2) {
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
  
  qX1 <- qr(X1)
  ## Take care of aliased variables and pivoting in rhs
  X1 <- X1[, qX1$pivot, drop=FALSE]
  X1 <- X1[, 1:qX1$rank, drop=FALSE]
  
  p <- NCOL(X1)
  mod <- lm(y ~ X1-1)
  s2  = summary(mod)$s**2
  D0  = diag(n)#diag(mu*(1-mu)) for logistic
  res = resid(mod)
  P0  = D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)
  px  = NCOL(X1)
  
  if (is.matrix(Ks)){Ks = list(Ks)}
  Qs = lapply(Ks, getQ, res, s2)  
  if (method == "davies"){    
    lambda0 = lapply(Ks, getLambda_davies, P0)
    if (length(Ks) == 1){
      p = getindivP_davies(Qs[[1]], lambda0[[1]], n, px)    
      return(list(indivP = p, KRV = KRVs, R2=R2))  # This should prevent the algorithm going any further
    }
    ps = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      ps[i] = getindivP_davies(Qs[[i]], lambda0[[i]], n, px)    
    }
    minP    = min(ps)
    
    #     q_sim = sapply(1:length(Ks), function(j){
    #       sapply(1:nperm, function(i) {
    #         ind <- sample(n)
    #         Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
    #       })  
    #     }) 
    
    q_sim =  sapply(1:nperm, function(i) {
      ind <- sample(n)
      sapply(1:length(Ks),function(j){
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
      })
    })  
    
    q_sim = t(q_sim)
    
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1) 
    p_perm = p_all[1,]
    
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  
    
    return(list(indivP = ps , omnibus_p = p_final, KRV = KRVs, R2 = R2))    
  } 
  if (method == "moment"){ 
    parm = sapply(Ks, getParamSat, P0, px)  
    if (length(Ks) == 1){
      p_sat = getSat(Qs[[1]], keppa_tlt = parm[,1]$keppa_tlt, niu_tlt = parm[,1]$niu_tlt) 
      return(list(indivP = as.numeric(unlist(p_sat)), KRV = KRVs, R2 = R2))
    }
    p_sat = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      p_sat[i] = getSat(Qs[[i]], keppa_tlt = parm[,i]$keppa_tlt, niu_tlt = parm[,i]$niu_tlt)    
    }
    minP = min(p_sat)
    #     
    #     q_sim = sapply(1:length(Ks), function(j){
    #       sapply(1:nperm, function(i) {
    #         ind <- sample(n)
    #         Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) # the adjusted is zero in this case 
    #       })  
    #     }) 
    #     
    
    q_sim =  sapply(1:nperm, function(i) {
      ind <- sample(n)
      sapply(1:length(Ks),function(j){
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
      })
    })  
    q_sim = t(q_sim)
    
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank) -1)/(nperm +1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
    p_perm = p_all[1,]
    
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm + 1)
    
    return(list(indivP = p_sat, omnibus_p = p_final, KRV = KRVs, R2 = R2)) 
  }
  if (method == "permutation"){
    #     
    #     q_sim = sapply(1:length(Ks), function(j){
    #       sapply(1:nperm, function(i) {
    #         ind <- sample(n)
    #         Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) # the adjusted is zero in this case 
    #       })  
    #     }) 
    
    
    q_sim =  sapply(1:nperm, function(i) {
      ind <- sample(n)
      sapply(1:length(Ks),function(j){
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
      })
    })  
    
    
    if (length(Ks) == 1){
      p_perm = (sum(q_sim > Qs)+ 1)/(nperm + 1) 
      return(list(indivP = p_perm, KRV = KRVs, R2 = R2))
      
    }
    q_sim = t(q_sim)
    
    # P = (# of statistics permutation >= observed test statistics + 1) / (#permutation + 1), 
    # thus the smallest p value will be  1/ (# permutation + 1).
    
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank) -1)/(nperm +1)  
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm + 1)
    
    
    # = (# of statistics permutation >= observed test statistics + 1) / (#permutation + 1)
    
    return(list(indivP = p_perm , omnibus_p = p_final, KRV = KRVs, R2 = R2)) 
    
  }
  
}
