MiRKATS <-
function(kd, distance=FALSE, obstime, delta, covar=NULL, beta=NULL, perm=FALSE, nperm=1000){
  # light checks for input
  if(length(obstime) != length(delta)) stop("Please make sure you have n observed times and n event indicators.")
  if(!is.null(beta) & is.null(covar)) warning("Your input includes coefficients but no covariates. Did you intend to include covariates?")
  if(nrow(kd) != length(obstime)) stop("Number of observed times does not match distance or kernel matrix. Please check your object dimensions.")
  if(class(kd) != "matrix") stop("Please convert your kernel or distance object into a matrix.")
  if(length(obstime) >= 100 & perm==TRUE){warning("Permutation p-values are not recommended unless n<100. Computation time may be long.")}
  if(length(obstime) <= 50 & perm==FALSE){warning("May wish to use permutation p-values when n <= 50.")}

  # sort in order of observed times
  ord = order(obstime)
  n = length(obstime)
  U = obstime[ord]
  D = delta[ord]
  if(!is.null(covar)){ X = as.matrix(as.matrix(covar, nrow=n)[ord,], nrow=n) } else { X = NULL }
  if(distance==TRUE){ kk <- D2K(kd); K <- kk[ord,ord] }
  if(distance==FALSE){ K = kd[ord,ord] }
  ftimes <- sort(U[D==1])     # failure times (i.e., event observed)
  fdups <- table(ftimes)
  fdups <- unlist(sapply(fdups, FUN=function(x)rep(x, x))) # counts at each failure time
  fdupcts <- unlist(sapply( table(ftimes), FUN=function(x) c(0:(x-1)) ))

                                        # if no covariates or no coefficients provided
  if(!is.null(X) & is.null(beta)){ beta <- coxph(Surv(U,D) ~ ., data=data.frame(X), ties="efron")$coef }
  if(is.null(X)){ eta <- rep(0,n) } else{ eta = X %*% beta }
  V <- as.vector(exp(eta))

  # notation (sometimes) matches Lin 2011
  # Efron approximation for ties
  R.i <- sapply(ftimes, FUN=function(x) U >= x)  ## risk sets (n x # failure times)
  D.i <- sapply(ftimes, FUN=function(x) U==x & D==1)  ## events at time t_i
  num <- apply(R.i, 2, FUN=function(x) V*x) - t(apply( apply(D.i, 2, FUN=function(x) x*V), 1, FUN=function(x) x*(fdupcts/fdups) ))
  denom <- apply( apply(R.i, 2, FUN=function(x) V*x) - t(apply( apply(D.i, 2, FUN=function(x) x*V), 1, FUN=function(x) x*(fdupcts/fdups) )), 2, sum)
  dLam <- apply(num, 1, FUN=function(x) x/denom)
  Lam <- apply(dLam, 2, sum)
  M.hat <- D - Lam  ## matches Martingale residuals from coxph

  if(perm==TRUE){
    obs.stat <- as.numeric(M.hat%*%K%*%M.hat)
    stat <- c()
    count <- 0

    while(count < nperm){
      # permute
      ord <- sample(n)
      mhat <- M.hat[ord]
      this <- mhat%*%K%*%mhat

      stat <- c(stat, this)
      count <- count + 1
    }
    pval <- mean(stat >= obs.stat)
  } else {
    # P0.star (from Han Chen)
    V.mat <- diag(Lam) - t(dLam) %*% dLam
    if(is.null(X)){ P0 <- V.mat } else { P0 <- V.mat - V.mat%*%X %*% solve(t(X)%*%V.mat%*%X) %*% t(X)%*%V.mat }
    eig.P0 <- eigen(P0)
    eig.P0.val <- eig.P0$values
    eig.P0.val[eig.P0.val < 1e-12] <- 0   # round eigenvalues too small to take the sqrt; enforce psd
    sqrt.P0 <- eig.P0$vectors %*% diag(sqrt(eig.P0.val)) %*% t(eig.P0$vectors)

    # Small sample correction
    dLam.2 <- dLam^2
    Lam.2 <- apply(dLam.2, 2, sum)
    W <- Lam - Lam.2
    if(sum(which(is.na(W)))>0) { W[which(is.na(W))] <- 1e-10 }
    W[which(W<1e-12)] <- 1e-10

    W.mat <- diag(W)
    this.eig <- eigen(W.mat)
    sqrt.W <- this.eig$vectors %*% diag(sqrt(this.eig$values)) %*% t(this.eig$vectors)
    z <- eta + (1/W)*M.hat
    if(!is.null(X)){
      y.star <- sqrt.W %*% z
      X.star <- sqrt.W %*% X
      e.star <- y.star - X.star %*% beta
      P0.star <- diag(rep(1,n)) - X.star %*% solve(t(X.star) %*% X.star) %*% t(X.star)
    } else {
      P0.star <- diag(n)
    }

    # Calculate test statistic and p-value
    obs.r <- c( M.hat %*% K %*% M.hat / (sum(M.hat^2)) )
    stat1 <- P0.star %*% (K - obs.r*diag(n)) %*% P0.star
    stat <- sqrt.P0 %*% stat1 %*% sqrt.P0
    lam <- eigen(stat)$values
    pval <- davies(0, lam)$Qq
  }
  return(pval)
}
