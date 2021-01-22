#' Small-sample SKAT for correlated (continuous) data ('c' stands for 'correlated'). Called within GLMM-MiRKAT. 
#'
#' Compute the adjusted score statistic and p-value
#' @param lmer.obj  A fitted lme4 object (model under H0)
#' @param Ks    A kernel matrix or list of kernels, quantifying the similarities between samples. 
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
CSKAT <- function(lmer.obj, Ks){
  
  no.omnibus <- FALSE 
  
  # Prepare data 
  if (is.matrix(Ks)) {
    Ks = list(Ks) 
    no.omnibus = TRUE
  } else if (!is.list(Ks)) {
    stop("Please enter either a single kernel matrix or a list of kernels for Ks.")
  } 
  
  # Run CSKAT for all candidate kernels using Davies for p-value calculation 
  pvs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    # returns q-statistic, p-value, variance matrix, residuals 
    res <- inner.CSKAT(lmer.obj, K = Ks[[j]]) 
    pvs[j] <- res$p.value 
  }
  names(pvs) = names(Ks)
  
  
  if(!no.omnibus){
    cauchy.t <- sum(tan((0.5 - pvs)*pi))/length(pvs)
    pv.opt <- 1 - pcauchy(cauchy.t)
  }
    
  if (!no.omnibus) {
    return(list(p_values = pvs, omnibus_p = pv.opt)) 
  } else { 
    return(p_values = pvs) 
  }
} 
