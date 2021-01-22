#'The Microbiome Regression-based Kernel Association Test Based on the Generalized Linear Mixed Model
#' 
#' GLMMMiRKAT utilizes a generalized linear mixed model to allow dependence among samples. 
#' 
#' Missing data is not permitted. Please remove all individuals with missing y, X, and Ks prior to input for analysis.
#' 
#' y and X (if not NULL) should be numerical matrices or vectors with the same number of rows. 
#'
#' Ks should either be a list of n by n kernel matrices (where n is sample size) or a single kernel matrix. If you have distance
#' matrices from metagenomic data, each kernel can be constructed through function D2K. Each kernel can also be constructed 
#' through other mathematical approaches.
#'          
#' If model="gaussian" and method="davies", CSKAT is called. CSKAT utilizes the same omnibus test as GLMMMiRKAT. See ?CSKAT for more details.
#' 
#' The "method" argument only determines kernel-specific p-values are generated. When Ks is a list of multiple kernels,
#' an omnibus p-value is computed via permutation.
#' 
#' 
#' @param y A numeric vector of Gaussian (e.g., body mass index), Binomial (e.g., disease status, treatment/placebo) or 
#' Poisson (e.g., number of tumors/treatments) traits.
#' @param X A vector or matrix of numeric covariates, if applicable (default = NULL).
#' @param Ks A list of n-by-n OTU kernel matrices or one singular n-by-n OTU kernel matrix, where n is sample size.
#' @param id A vector of cluster (e.g., family or subject including repeated measurements) IDs. Defaults to NULL since it is 
#' unnecessary for the CSKAT call.
#' @param time.pt A vector of time points for the longitudinal studies. 'time.pt' is not required (i.e., 'time.pt = NULL') 
#' for the random intercept model. Default is time.pt = NULL. 
#' @param model A string declaring which model ("gaussian", "binomial" or "poisson") is to be used; should align with whether a 
#' Gaussian, Binomial, or Poisson trait is being inputted for the y argument.
#' @param method A string declaring which method ("perm" or "davies) will be used to calculate the p-value. Davies is only 
#'  available for Gaussian traits. Defaults to "perm". 
#' @param slope An indicator to include random slopes in the model (slope = TRUE) or not (slope = FALSE). 'slope = FALSE' is for 
#' the random intercept model. 'slope = TRUE' is for the random slope model. For the random slope model (slope = TRUE), 'time.pt' 
#' is required.
#' @param omnibus A string equal to either "Cauchy" or "permutation" (or nonambiguous abbreviations thereof), specifying whether 
#'  to use the Cauchy combination test or residual permutation to generate the omnibus p-value. 
#' @param nperm The number of permutations used to calculate the p-values and omnibus p-value. Defaults to 5000.
#' 
#' @return Returns a p-value for each inputted kernel matrix, as well as an overall omnibus p-value if more than one kernel matrix
#'         is inputted
#'         \item{p_values}{p-value for each individual kernel matrix}
#'         \item{omnibus_p}{overall omnibus p-value calculated by permutation for the adaptive GLMMMiRKAT analysis}
#'         
#'@author 
#'Hyunwook Koh
#'
#'@references
#'Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed 
#'model for correlated microbiome studies. Front. Genet. 458(10), 1-14.
#'
#'
#'@importFrom lme4 lmer glmer getME 
#'@importFrom Matrix Diagonal 
#'@import permute 
#'
#'
#'@examples
#'
#'library(vegan) 
#'
#'## Example with Gaussian (e.g., body mass index) traits
#'## For non-Gaussian traits, see vignette. 
#'
#'# Import example microbiome data with Gaussian traits
#'data(nordata)
#'otu.tab <- nordata$nor.otu.tab
#'meta <- nordata$nor.meta
#'
#'# Create kernel matrices
#'# could use phylogenetic kernels as below; computation time is slightly higher
#'# tree <- nordata$nor.tree
#'# unifracs <- GUniFrac::GUniFrac(otu.tab, tree, alpha=c(1))$unifracs
#'D_BC = as.matrix(vegdist(otu.tab, 'bray'))
#'K_BC = D2K(D_BC)
#'
#'# Run GLMM-MiRKAT
#'GLMMMiRKAT(y = meta$y, X = cbind(meta$x1, meta$x2), id = meta$id, 
#'           Ks = K_BC, model = "gaussian", nperm = 500)
#'
#' @export
GLMMMiRKAT <- function(y, X = NULL, Ks, id = NULL, time.pt = NULL, model, method = "perm" , 
                       slope = FALSE, omnibus = "perm", nperm = 5000){
  
  om <- substring(tolower(omnibus), 1, 1)
  
  if (is.null(time.pt) & slope) {
    stop("time.pt is required for the random slope model")
  }
  if(model != "gaussian" & method == "davies"){
    stop("Davies is only available for Gaussian traits")
  }
  
  if (model == "gaussian" & method =="davies"){  ## Calls CSKAT 
    if (om == "p") {
      warning("Permutation omnibus test not available with Davies P-values. Defaulting to Cauchy combination test.")
    }
    y <- scale(y)
    if (is.null(X)) {
      if (!is.null(time.pt) & slope) {
        fit <- lmer(y ~ (time.pt | id))
      } else {
        fit <- lmer(y ~ (1 | id))
      }
    } else {  # if X not null 
      if (!is.null(time.pt) & slope) {
        fit <- lmer(y ~ (time.pt | id) + X)
      } else {
        fit <- lmer(y ~ (1 | id) + X)
      }
    }
    out <- CSKAT(fit, Ks=Ks)
    return(out)  #This should stop the function from going any further 
  }
  
  if (is.matrix(Ks)) {
    Ks <- list(Ks) 
    no.omnibus <- TRUE 
  } else {
    no.omnibus <- FALSE 
  }
  
  if (is.null(time.pt)) {
    id <- as.character(id)
    ind <- order(id)
    id <- id[ind]
    y <- y[ind]
    if (is.matrix(X)|is.data.frame(X)) {
      X <- as.matrix(X)[ind, ] 
    } else {
      X <- X[ind]
    }
    
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
  } else {
    ind <- order(as.character(time.pt))
    id <- as.character(id)
    id <- id[ind]
    y <- y[ind]
    if (is.matrix(X)|is.data.frame(X)) X <- as.matrix(X)[ind, ] else X <- X[ind]
    time.pt <- time.pt[ind]
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
    ind <- order(id)
    id <- id[ind]
    y <- y[ind]
    if (is.matrix(X)|is.data.frame(X)) X <- as.matrix(X)[ind, ] else X <- X[ind]
    time.pt <- time.pt[ind]
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
  }
  
  if (model == "gaussian") {
    y <- scale(y)
    if (is.null(X)) {
      if (!is.null(time.pt) & slope) {
        fit <- lmer(y ~ (time.pt | id))
      } else {
        fit <- lmer(y ~ (1 | id))
      }
    } else {  # if X not null 
      if (!is.null(time.pt) & slope) {
        fit <- lmer(y ~ (time.pt | id) + X)
      } else {
        fit <- lmer(y ~ (1 | id) + X)
      }
    }
  } else {  # if not Gaussian 
    if (is.null(X)) {
      if (!is.null(time.pt) & slope) {
        fit <- glmer(y ~ (time.pt | id), family = model)
      } else {
        fit <- glmer(y ~ (1 | id), family = model)
      }
    } else {
      if (!is.null(time.pt) & slope) {
        fit <- glmer(y ~ (time.pt | id) + X, family = model)
      } else {
        fit <- glmer(y ~ (1 | id) + X, family = model)
      }
    }
  }
  
  fe <- as.numeric(getME(fit, "X") %*% getME(fit, "fixef"))
  re <- as.numeric(getME(fit, "Z") %*% getME(fit, "b"))
  if (model == "gaussian") {
    mu <- fe + re
    inv.de.mu <- diag(length(y))
  }
  if (model == "binomial") {
    mu <- exp(fe + re)/(1 + exp(fe + re))
    inv.de.mu <- diag(as.numeric(mu * (1 - mu)))
  }
  if (model == "poisson") {
    mu <- exp(fe + re)
    inv.de.mu <- diag(as.numeric(mu))
  }
  y.star <- as.numeric(fe + re + ginv(inv.de.mu) %*% (y - mu))
  r <- y.star - fe
  r.X <- as.matrix((getME(fit, "Z") %*% Matrix::crossprod(getME(fit, "Lambdat")) %*% getME(fit, "Zt")))
  if (model == "gaussian") {
    e.var <- diag(rep(var(mu), length(y)))
  }
  if (model == "binomial") {
    e.var <- diag(rep(1, length(y)))
  }
  if (model == "poisson") {
    e.var <- diag(rep(1, length(y)))
  }
  v.X <- r.X + e.var
  inv.vX <- ginv(v.X)
  
  r.s <- list()
  if (!slope) {
    id.time.pt <- id
    p.ind <- as.matrix(getPermuteMatrix(nperm, length(r), 
                                        strata = id.time.pt))
    clust.sizes <- as.numeric(names(table(table(id.time.pt))))
    block.r.ids <- list()
    for (l in 1:nperm) {
      block.r.id <- list()
      for (j in 1:length(clust.sizes)) {
        block.r.id.clust <- list()
        u.id.clust <- names(which(table(id.time.pt) == clust.sizes[j]))
        r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
        for (k in 1:length(u.id.clust)) {
          block.r.id.clust[[k]] <- which(id.time.pt == r.u.id.clust[k])
        }
        block.r.id[[j]] <- unlist(block.r.id.clust)
      }
      unlist.block.r.id <- rep(NA, length(id.time.pt))
      for (j in 1:length(clust.sizes)) {
        unlist.block.r.id[id.time.pt %in% names(which(table(id.time.pt) == 
                                                        clust.sizes[j]))] <- block.r.id[[j]]
      }
      block.r.ids[[l]] <- unlist.block.r.id
    }
    for (l in 1:nperm) {
      p.ind[l, ] <- p.ind[l, block.r.ids[[l]]]
    }
    for (j in 1:nperm) {
      r.s[[j]] <- r[p.ind[j, ]]
    }
  }
  if (!is.null(time.pt) & slope) {
    id.names <- names(table(id))
    id.time.pt <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      id.time.pt[ind] <- paste(c(id.names[j], as.character(time.pt[ind])), 
                               collapse = ".")
    }
    noid.time.pt <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      noid.time.pt[ind] <- paste(as.character(time.pt[ind]), 
                                 collapse = ".")
    }
    ex.clust <- names(table(noid.time.pt))
    
    ids <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      ids[ind] <- id.names[j]
    }
    if (length(names(which(table(noid.time.pt)==1))) != 0) {
      warn.ids <- unique(ids[noid.time.pt %in% names(which(table(noid.time.pt)==1))])
    }
    if (length(warn.ids) == 1) warning("The cluster id (", warn.ids, ") is a silent block which is not exchangeable with any other blocks. It did not contribute to the analysis.")
    if (length(warn.ids) > 1) warning("The cluster ids (", paste(warn.ids, collapse=", "), ") are silent blocks which are not exchangeable with any other blocks. They did not contribute to the analysis.")
    block.r.ids <- list()
    for (l in 1:nperm) {
      block.r.id <- list()
      for (j in 1:length(ex.clust)) {
        block.r.id.clust <- list()
        u.id.clust <- unique(id.time.pt[which(noid.time.pt == ex.clust[j])])
        r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
        for (k in 1:length(u.id.clust)) {
          block.r.id.clust[[k]] <- which(id.time.pt == r.u.id.clust[k])
        }
        block.r.id[[j]] <- unlist(block.r.id.clust)
      }
      unlist.block.r.id <- rep(NA, length(id.time.pt))
      for (j in 1:length(ex.clust)) {
        unlist.block.r.id[id.time.pt %in% unique(id.time.pt[which(noid.time.pt == ex.clust[j])])] <- block.r.id[[j]]
      }
      block.r.ids[[l]] <- unlist.block.r.id	
    }
    for (j in 1:nperm) {
      r.s[[j]] <- r[block.r.ids[[j]]]
    }
  }
  
  Qs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    Qs[j] <- (t(r) %*% inv.vX %*% Ks[[j]] %*% inv.vX %*% 
                r)/(t(r) %*% inv.vX %*% r)
  }
  Q0s <- list()
  for (j in 1:length(Ks)) {
    Q0s.inv <- rep(NA, nperm)
    for (k in 1:nperm) {
      Q0s.inv[k] <- (t(r.s[[k]]) %*% inv.vX %*% Ks[[j]] %*% 
                       inv.vX %*% r.s[[k]])/(t(r.s[[k]]) %*% inv.vX %*% 
                                               r.s[[k]])
    }
    Q0s[[j]] <- Q0s.inv
  }
  pvs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    pvs[j] <- (length(which(Q0s[[j]] > Qs[[j]])) + 1)/(nperm + 1)
  }
  
  if (!no.omnibus) {
    if (om == "p") {
      T <- min(pvs)
      T0 <- rep(NA, nperm)
      for (l in 1:nperm) {
        T0.s.n <- list()
        for (m in 1:length(Ks)) {
          T0.s.n[[m]] <- Q0s[[m]][-l]
        }
        a.Ts <- unlist(lapply(Ks, function(x) return((t(r.s[[l]]) %*% 
                                                        inv.vX %*% x %*% inv.vX %*% r.s[[l]])/(t(r.s[[l]]) %*% 
                                                                                                 inv.vX %*% r.s[[l]]))))
        a.pvs <- unlist(mapply(function(x, y) (length(which(x > 
                                                              y)) + 1)/nperm, T0.s.n, a.Ts))
        T0[l] <- min(a.pvs)
      }
      pv.opt <- (length(which(T0 < T)) + 1)/(nperm + 1)
    } else if (om == "c") {
      cauchy.t <- sum(tan((0.5 - pvs)*pi))/length(pvs)
      pv.opt <- 1 - pcauchy(cauchy.t)
    } else {
      stop("I don't know that omnibus option. Please choose 'permutation' or 'Cauchy'.")
    }
  }
  
  names(pvs) <- names(Ks)
  
  if (no.omnibus) {
    return(list(p_values = pvs))
  } 
  return(list(p_values = pvs, omnibus_p = pv.opt))
}

