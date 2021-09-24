# This is the file according to quantreg, but adjusted for LiQRAT.
# It takes the sparse matrix dsol[,1:(J-1)] - dsol[,2:J]


NewRanks <- function (v, score = "wilcoxon", Delta = NULL)
{
  A2 <- 1
  taus = v$taus
  J = v$J
  dt = taus[2:J] - taus[1:(J-1)]
  if (score == "wilcoxon") {
    A2 <- 1/12
    phi <- c(0,  0.5*taus[2:(J-1)] - 0.5*(taus[2:(J-1)])^2 , 0)
    dphi <- diff(phi)
  }
  else if (score == "normal") {
    phi <- c(0, dnorm(qnorm(taus[2:(J-1)])), 0)
    dphi <- diff(phi)
  }
  else if (score == "normalscale") {
    A2 = 2
    Qt <- qnorm(taus[2:(J - 1)])
    phi <- c(0, -dnorm(Qt) * Qt, 0)
    dphi <- diff(phi)
  }
  else if (score == "halfnormalscale") {
    A2 = 1.25
    Qt <- qnorm(taus[2:(J - 1)])
    phi <- c(0, (taus[2:(J - 1)] - dnorm(Qt) * Qt) * (taus[2:(J -  1)] > 0.5), 1)
    dphi <- diff(phi)
    dphi[dphi > 0.5] <- 0
    
  }
  else if (score == "lehmann") {
    taus <- taus[2:(J - 1)]
    phi <- c(0, (1-taus) * log(1 - taus), 0)
    dphi <- diff(phi)
    
  }
  else if(score == "trimmedlehmann"){
    # taus <- taus[2:(J - 1)]
    ind1 = length(which(taus <= Delta))
    ind2 = length(which(taus > Delta))
    phi1 = -(taus - 1) * log(1 - taus)
    phi2 = 1/(1-Delta)*(0.5*taus^2)-taus*(log(1-Delta)+1+Delta/(1-Delta))
    dphi1 = diff(phi1[1:(ind1+1)])
    dphi2 = diff(phi2[(ind1+1):J])
    dphi <- c(dphi1, dphi2)
    A2 = (1-Delta)*log(1-Delta)+Delta +(1-Delta)/3 -(1-Delta)^2/4
    
  }
  else if(score == "inverselehmann"){
    #taus = 1-taus
    #dt = taus[2:J]- taus[1:(J-1)]
    taus <- taus[2:(J - 1)]
    phi <- c(0, taus * log(taus), 0)
    dphi <- diff(phi)
    
  }
  else if(score == "trimmedinverselehmann"){
    # taus <- taus[2:(J - 1)]
    ind1 = length(which(taus >= Delta))
    ind2 = length(which(taus < Delta))
    phi1 = taus * log(taus)
    phi2 = 1/(Delta)*(0.5*(taus)^2)+taus*(log(Delta))
    dphi1 = diff(phi1[(ind2-1):J])
    dphi2 = diff(phi2[1:(ind2-1)])
    dphi <- c(dphi2, dphi1)
    A2 = (Delta)*log(Delta)+(1-Delta) +(Delta)/3 -(Delta)^2/4
    
  }
  
  else stop("invalid score function")
  
  if(score == "trimmedlehmann"){
    ranks <- as.vector((v$diff_dsol) %*% (dphi/dt)) - (1-Delta)/2
  }else if(score == "trimmedinverselehmann"){
    ranks <- as.vector((v$diff_dsol) %*% (dphi/dt)) + (Delta)/2
  }else{
    ranks <- as.vector((v$diff_dsol) %*% (dphi/dt))
  }
  
  return(list(ranks = ranks, A2 = A2))
  
}


# Transfer quantreg v

transfer_v <- function(v, diff_dsol = FALSE){
  v.new <- list()
  if(diff_dsol){
    v.new$J = J <- ncol(v$beta)
    v.new$diff_dsol = v$diff_dsol
    v.new$taus = taus <- v$beta[1, ]
  }else{
    v.new$J = ncol(v$sol)
    v.new$diff_dsol = (v$dsol[, 2:v.new$J] - v$dsol[, 1:(v.new$J - 1)])
    v.new$taus = v$sol[1, ]
  }
  return(v = v.new)
}
