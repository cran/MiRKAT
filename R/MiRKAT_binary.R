#' Microbiome Regression-Based Kernel Association Test for binary outcomes
#' 
#' Called by MiRKAT if the outcome variable is dichotomous (out_type="D")
#'  
#' This function is called by the exported function MiRKAT if the argument "out_type" of MiRKAT is equal to "D" (for dichotomous).
#'  
#' Each argument of MiRKAT_continuous is given the value of the corresponding argument given by the user to MiRKAT.
#'  
#' Function not exported
#'  
#' @param y A numeric vector of the dichotomous outcome variable 
#' @param X A numerical matrix or data frame, containing additional covariates that you want to adjust for (Default = NULL). If it is 
#' NULL, a intercept only model was fit.
#' @param Ks A list of n by n kernel matrices (or a single n by n kernel matrix), where n is the sample size. It can be constructed 
#' from microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.
#' @param method A string telling R which method to use to compute the kernel specific p-value (default = "davies"). "davies" 
#' represents an exact method that computes the p-value by inverting the characteristic function of the mixture chisq. We adopt an 
#' exact variance component tests because most of the studies concerning microbiome compositions have modest sample size. "moment" 
#' represents an approximation method that matches the first two moments. "permutation" represents a permutation approach for 
#' p-value calculation. 
#' @param family A string describing the error distribution and link function to be used in the linear model.
#' @param nperm the number of permutations if method = "permutation" or when multiple kernels are considered. if method = "davies" 
#'  or "moment", nperm is ignored.
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#'
#'  
#' @return If only one candidate kernel matrix is considered, returns a list containing the p-value for the candidate kernel matrix. 
#' If more than one candidate kernel matrix is considered, returns a list with two elements: the individual p-values for
#' each candidate kernel matrix, and the omnibus p-value. 
#'     \item{indivP}{p-value for each candidate kernel matrix}
#'     \item{omnibus_p}{omnibus p-value if multiple kernel matrices are considered}
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'
#'
#'
#'@author Ni Zhao
#'
#'@references 
#' Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J., Ringel, Y., Li, H. and Wu, M.C. (2015)).
#'  Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807
#'
#' Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and 
#' Microbiome Association Studies. 40: 5-19. doi: 10.1002/gepi.21934
#'
#' Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal
#'  Statistical Society. Series C , 29, 323-333.
#'
#' Satterthwaite, F. (1946). An approximate distribution of estimates of variance components. Biom. Bull. 2, 110-114.
#'
#' Lee S, Emond MJ, Bamshad MJ, Barnes KC, Rieder MJ, Nickerson DA; NHLBI GO Exome Sequencing Project-ESP Lung Project Team, 
#' Christiani DC, Wurfel MM, Lin X. (2012) Optimal unified approach for rare variant association testing with application to 
#' small sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.
#'
#' Zhou, J. J. and Zhou, H.(2015) Powerful Exact Variance Component Tests for the Small Sample Next Generation Sequencing Studies 
#' (eVCTest), in submission.



MiRKAT_binary = function(y, X = NULL, Ks, method = "davies", family= "binomial", nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  
  n <- length(y)
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
  
  if (returnKRV | returnR2) {
    reskrv = scale(resid(glm(y ~ 1, family = "binomial")))
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
  options(warn=2)  # make sure this model is correct
  mod <- glm(y ~ X1-1, family = family)
  options(warn=1)
  
  px  = NCOL(X1)
  mu  = mod$fitted.values
  res = y - mu  
  
  w   = mu*(1-mu)
  D0  =  sqrt(w)  
  DX12 = D0 * X1
  P0 = diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)
  if (method == "davies"){    
    if (n < 50){
      warning("For binary outcome and n < 50, p-value using davies method can be inaccurate at tails, permutation is recommended.")
    }
    S = sapply(Ks, getIndivP_binary, res,  D0, px, P0)
    ps = as.numeric(unlist(S[3,]))
    if (length(Ks) ==1){
      return(list(indivP = ps, KRV = KRVs, R2 = R2))
    }
    eP0 = c(rep(1, n-px), rep(0, px))
    Qs = unlist(S[1,])
    
    # Is it so difficult to do that permutation in this case? 
    
    #     p_sim = sapply(1:length(Ks), function(j){
    #       p1 = sapply(1:nperm, function(i) {
    #         ind <- sample(n)
    #         Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
    #         PKP = P0 %*% (D0*t(D0 * Ks[[j]][ind, ind])) %*% P0
    #         ePKP = eigen(PKP, symmetric = T)$values
    #         lambda0 = ePKP - Q1*eP0/n
    #         lambda0 = lambda0[abs(lambda0) >= 1e-10]
    #         p = davies(0, lambda=lambda0, acc = 1e-6)$Qq
    #         p = ifelse(p < 0, 0, p)
    #         return(p)
    #       })  
    #     }) 
    
    
    q_sim = sapply(1:nperm, function(i){
      ind <- sample(n)
      p1 = sapply(1:length(Ks), function(j) {
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
        return(Q1)
      })  
    }) 
    
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  
    
    # Naming the p-values with the names of their corresponding kernel matrices
    out_pvs = ps 
    names(out_pvs) = names(Ks) 
  } else if (method == "moment"){
    S = sapply(Ks, getIndivP_hm, res, mu, D0, P0)
    ps = as.numeric(unlist(S[1,]))
    Qs = unlist(S[2,])
    
    if (length(Ks) == 1){
      return(list(indivP = ps, KRV = KRVs, R2 = R2))
    }
    q_sim = sapply(1:nperm, function(i){
      ind <- sample(n)
      p1 = sapply(1:length(Ks), function(j) {
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
        return(Q1)
      })  
    }) 
    
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  
    
    out_pvs = ps 
    names(out_pvs) = names(Ks) 
  } else if (method == "permutation"){
    
    Qs = lapply(Ks, getQ, res, s2 = 1)
    
    q_sim = sapply(1:nperm, function(i){
      ind <- sample(n)
      p1 = sapply(1:length(Ks), function(j) {
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
        return(Q1)
      })  
    }) 
    
    if (length(Ks) == 1){
      p_perm = (sum(q_sim > Qs)+ 1)/(nperm + 1) 
      #Naming the outputted p-values with the names of their corresponding kernel matrices
      names_Ks <- names(Ks)
      named_pvals <- setNames(p_perm, names_Ks)
      return(list(indivP = named_pvals, KRV = KRVs, R2 = R2))
    }
    
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)   
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm + 1)
    
    # Naming the outputted p_values with the names of their correspinding kernel matrices
    out_pvs = p_perm 
    names(out_pvs) = names(Ks)
  }
  return(list(p_values = out_pvs, omnibus_p = p_final, KRV = KRVs, R2 = R2)) 
}
