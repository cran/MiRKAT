#' Small-sample SKAT for correlated (continuous) data ('c' stands for 'correlated') 
#'
#' Compute the adjusted score statistic and p-value
#' @param formula.H0  A two-sided linear formula object describing both the fixed-effects and random-effects part of
#'  the model under the null, use the same syntax as the "lmer" in "lme4" package 
#' @param data  An optional data frame containing the variables named in formula. Default: NULL. 
#' @param Ks    A kernel matrix or list of kernels, quantifying the similarities between samples. 
#' @param nperm  Number of permutations for calculating the omnibus p-value. Ignored unless Ks is a list of candidate kernels. 
#' 
#' @return
#' \describe{
#'   \item{p.value}{Association p-values}
#'   \item{Q.adj}{Adjusted score statistics}
#' }
#' @author 
#' Nehemiah Wilson, Anna Plantinga, Xiang Zhan, Jun Chen. 
#' 
#' @references
#' Zhan X, et al. (2018) A small-sample kernel association test for correlated data with application to microbiome association studies.  
#' Genet Epidemiol.
#' 
#' @examples
#' 
#' Y <- rnorm(100)
#' Z <- matrix(rnorm(200), 100, 2)
#' ID <- gl(20, 5)
#' G <- matrix(rbinom(1000, 2, 0.05), 100, 10)
#' K <- G %*% t(G)
#' CSKAT(formula.H0 = Y ~ Z + (1 | ID), Ks = K)
#' 
#' @export 
#' 
CSKAT <- function(formula.H0, data = NULL, Ks, nperm = 999){
  
  # Prepare data 
  if (is.matrix(Ks)) {
    Ks = list(Ks) 
  } else if (!is.list(Ks)) {
    stop("Please enter either a single kernel matrix or a list of kernels for Ks.")
  } 
  
  # Run CSKAT for all candidate kernels using Davies for p-value calculation 
  pvs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    # returns q-statistic, p-value, variance matrix, residuals 
    res <- inner.CSKAT(formula.H0 = formula.H0, data = data, K = Ks[[j]]) 
    pvs[j] <- res$p.value 
  }
  # same resids under H0 regardless of kernel 
  resid <- res$resid
  
  
  if(length(Ks) > 1){
    # GLMMMiRKAT Omnibus Test
    T <- min(pvs) # test statistic for GLMMMiRKAT omnibus test
    r.s <- list() # list containing nperm number of permuted residuals 
    for(i in 1:nperm){
      r.s[[i]] <- resid[shuffle(length(resid))]
    }
    
    Q0s <- list() 
    for (j in 1:length(Ks)) {
      s2 <- sum(resid^2)
      Q0s.inv <- rep(NA, nperm)
      for (k in 1:nperm) {
        Q0s.inv[k] <- (as.numeric(r.s[[k]] %*% Ks[[j]] %*% r.s[[k]] / s2 ))
      }
      Q0s[[j]] <- Q0s.inv
    }
    
    # Null p-values
    T0 <- rep(NA, nperm)
    for (l in 1:nperm) { 
      T0.s.n <- list()
      for (m in 1:length(Ks)) {
        T0.s.n[[m]] <- Q0s[[m]][-l]
      }
      
      
      a.Ts <- unlist(lapply(Ks, function(x) return(r.s[[l]] %*% x %*% r.s[[l]]))) #GLMMMiRKAT uses this...x is in the form of a list tho
      #so wont this not work? Must it be in a for loop to reach
      #each individual element in x?
      
      a.pvs <- unlist(mapply(function(x, y) return(length(which(x > y)) + 1)/nperm, T0.s.n, a.Ts))
      T0[l] <- min(a.pvs)
    }
    
    pv.opt <- (length(which(T0 < T)) + 1)/(nperm + 1)
  
    if(is.null(names(Ks))){
      warning("P-values are not labeled with their corresponding kernel matrix. In order to have them labeled,
              make your list of kernel matrices for the input of the form 'list(name1=K1, name2=K2'...) in order for the output
              p-values to be labeled with 'name1,' 'name2,' etc.")
    } else {
      names(pvs) = names(Ks)
    }
    
    return(list(p_values = pvs, omnibus_p = pv.opt)) 
  } else { 
    return(p_values = pvs) 
  }
} 
