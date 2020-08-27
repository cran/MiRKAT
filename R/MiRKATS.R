#'Microiome Regression-based Kernel Association Test for Survival
#' 
#' Community level test for association between microbiome composition and survival outcomes (right-censored time-to-event data) 
#' using kernel matrices to compare similarity between microbiome profiles with similarity in survival times. 
#' 
#' obstime, delta, and X should all have n rows, and the kernel or distance matrix should be a single n by n matrix. 
#' If a distance matrix is entered (distance=TRUE), a kernel matrix will be constructed from the distance matrix.
#' 
#' Update in v1.1.0: MiRKATS also utilizes the OMiRKATS omnibus test if more than one kernel matrix is provided by the user. 
#' The OMiRKATS omnibus test calculates an overall p-value for the test via permutation. 
#'
#' Missing data is not permitted. Please remove individuals with missing data on y, X or in the kernel or distance matrix prior 
#' to using the function.
#'
#' The Efron approximation is used for tied survival times.
#' 
#' @param obstime A numeric vector of follow-up (survival/censoring) times.
#' @param delta Event indicator: a vector of 0/1, where 1 indicates that the event was observed for a subject (so "obstime" is 
#' survival time), and 0 indicates that the subject was censored.
#' @param X A vector or matrix of numeric covariates, if applicable (default = NULL).
#' @param Ks A list of or a single numeric n by n kernel matrices or matrix (where n is the sample size).
#' @param beta A vector of coefficients associated with covariates. If beta is NULL and covariates are present, coxph is used to 
#' calculate coefficients (default = NULL).
#' @param perm Logical, indicating whether permutation should be used instead of analytic p-value calculation (default=FALSE). 
#' Not recommended for sample sizes of 100 or more.
#' @param omnibus A string equal to either "Cauchy" or "permutation" (or nonambiguous abbreviations thereof), specifying whether 
#'  to use the Cauchy combination test or residual permutation to generate the omnibus p-value. 
#' @param nperm Integer, number of permutations used to calculate p-value if perm==TRUE (default=1000) and to calculate omnibus p-value if omnibus=="permutation". 
#' @param returnKRV A logical indicating whether to return the KRV statistic. Defaults to FALSE. 
#' @param returnR2 A logical indicating whether to return the R-squared coefficient. Defaults to FALSE.  
#'
#'@return 
#'Return value depends on the number of kernel matrices inputted. If more than one kernel matrix is given, MiRKATS returns two
#'items; a vector of the labeled individual p-values for each kernel matrix, as well as an omnibus p-value from the Optimal-MiRKATS
#'omnibus test. If only one kernel matrix is given, then only its p-value will be given, as no omnibus test will be needed.
#'    \item{p_values}{individual p-values for each inputted kernel matrix}
#'    \item{omnibus_p}{overall omnibus p-value}
#'     \item{KRV}{A vector of kernel RV statistics (a measure of effect size), one for each candidate kernel matrix. Only returned if returnKRV = TRUE}
#'     \item{R2}{A vector of R-squared statistics, one for each candidate kernel matrix. Only returned if returnR2 = TRUE}
#'     
#' @author 
#' Nehemiah Wilson, Anna Plantinga
#' 
#' @references 
#' Plantinga, A., Zhan, X., Zhao, N., Chen, J., Jenq, R., and Wu, M.C. MiRKAT-S: a distance-based test of association between
#' microbiome composition and survival times. Microbiome, 2017:5-17. doi: 10.1186/s40168-017-0239-9
#'
#' Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J., Ringel, Y., Li, H. and Wu, M.C. (2015)).
#' Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807
#' 
#' Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and
#' Microbiome Association Studies. 40:5-19. doi: 10.1002/gepi.21934
#' 
#' Efron, B. (1977) "The efficiency of Cox's likelihood function for censored data." Journal of the American statistical 
#' Association 72(359):557-565.
#' 
#' Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal 
#' Statistical Society Series C, 29:323-333
#' 
#' @importFrom survival coxph Surv
#' 
#'@examples 
#'
#'###################################
#'# Generate data
#'library(GUniFrac)
#'
#'# Throat microbiome data
#'data(throat.tree)
#'data(throat.otu.tab)
#'
#'unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
#'Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], 
#'          BC= as.matrix(vegdist(throat.otu.tab, method="bray"))) 
#'Ks = lapply(Ds, FUN = function(d) D2K(d))
#'
#'# Covariates and outcomes
#'covar <- matrix(rnorm(120), nrow=60)
#'S <- rexp(60, 3)   # survival time 
#'C <- rexp(60, 1)   # censoring time 
#'D <- (S<=C)        # event indicator
#'U <- pmin(S, C)    # observed follow-up time
#'
#'MiRKATS(obstime = U, delta = D, X = covar, Ks = Ks, beta = NULL)
#'
#'
#' @export
MiRKATS <- function(obstime, delta, X = NULL, Ks,  beta = NULL, perm=FALSE, omnibus="permutation", nperm=999, returnKRV = FALSE, returnR2 = FALSE){
  
  om <- substring(tolower(omnibus), 1, 1)
  
  if(!is.list(Ks)){
    if (!is.matrix(Ks)) stop("Please convert your kernel into a matrix.")
    Ks <- list(Ks)
  }
  
  if(!length(Ks)==1){
    if(is.null(names(Ks))){
      message("Your p-values are not labeled with their corresponding kernel matrix. In order to have them labeled,
              make your list of kernel matrices for the input of the form 'list(name1=K1, name2=K2'...) in order for the output
              p-values to be labeled with 'name1,' 'name2,' etc.")
    }
  }
  
  # light checks for input
  if(length(obstime) != length(delta)) stop("Please make sure you have n observed times and n event indicators.")
  if(!is.null(beta) & is.null(X)) warning("Your input includes coefficients but no covariates. Did you intend to include covariates?")
  if(nrow(Ks[[1]]) != length(obstime)) stop("Number of observed times does not match distance or kernel matrix. Please check your object dimensions.")
  if(length(obstime) >= 100 & perm==TRUE){warning("Permutation p-values are not recommended unless n<100. Computation time may be long.")}
  if(length(obstime) <= 50 & perm==FALSE){warning("Permutation p-values are recommendeded when n <= 50.")}
  
  if (returnKRV | returnR2) {
    resids = scale(residuals(coxph(Surv(time = obstime, event = delta) ~ 1), type="martingale"))
    L = resids %*% t(resids) 
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
  
  
  # Calculate individual p-values
  pvals <- c()
  for(i in 1:length(Ks)){
    pvals[i] <- inner.MiRKATS(obstime=obstime, delta=delta, covar = X, K = Ks[[i]], beta=beta, perm=perm, nperm=nperm)
  }
  names(pvals) <- names(Ks)
  
    
  if(length(Ks) > 1){
    
    if (om == "p") {
      ######################################
      # Optimal-MiRKATS omnibus test begins
      ######################################
      
      if (is.null(X)) {X <- rep(1, length(obstime))}
      r <- coxph(Surv(obstime, delta) ~ ., data=as.data.frame(X))$residuals 
      r.s <- list() # becomes a list of permuted vectors of residuals
      for (j in 1:nperm) {
        r.s[[j]] <- r[shuffle(length(r))]
      }
      
      T0s.mirkats <- list()
      for (j in 1:length(Ks)) {
        T0s.mirkats.inv <- rep(NA, nperm)
        for (k in 1:nperm) {
          T0s.mirkats.inv[k] <- t(r.s[[k]])%*%Ks[[j]]%*%r.s[[k]]
        }
        T0s.mirkats[[j]] <- T0s.mirkats.inv
      }
      
      Q.mirkats <- min(pvals) # The test statistic for OMiRKATS, Q.mirkats, is the minimum of the p-values from MiRKATS
      Q0.mirkats <- rep(NA, nperm) 
      for (l in 1:nperm) { # Creating a list of omnibus test statistics (minimum p-values from null distributions of test statistics)
        Q0.mirkats.s.n <- list()
        for (m in 1:length(Ks)) {
          Q0.mirkats.s.n[[m]] <- T0s.mirkats[[m]][-l]
        }
        a.Qs.mirkats <- unlist(lapply(Ks,function(x) return(t(r.s[[l]])%*%x%*%r.s[[l]])))
        a.pvs <- unlist(mapply(function(x,y)length(which(abs(x) >= abs(y)))/(nperm-1),Q0.mirkats.s.n,a.Qs.mirkats))
        Q0.mirkats[l] <- min(a.pvs)
      } 
      p.omirkats <- length(which(Q0.mirkats <= Q.mirkats))/nperm # The omnibus p-value
      
    } else if (om == "c") {
      cauchy.t <- sum(tan((0.5 - pvals)*pi))/length(pvals)
      p.omirkats <- 1 - pcauchy(cauchy.t)
    } else {
      stop("I don't know that omnibus option. Please choose 'permutation' or 'Cauchy'.")
    }
    
    
    ## return if multiple kernels 
    if (is.null(KRVs) & is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = p.omirkats))
    } else if (is.null(KRVs) & !is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = p.omirkats, R2 = R2))
    } else if (!is.null(KRVs) & is.null(R2)) {
      return(list(p_values = pvals, omnibus_p = p.omirkats, KRV = KRVs))
    } else {
      return(list(p_values = pvals, omnibus_p = p.omirkats, KRV = KRVs, R2 = R2))    
    }

  }

  ## return if only one kernel 
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

