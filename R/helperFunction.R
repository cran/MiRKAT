
# For continuous outcome
permuted.index = function (n){
  out = sample.int(n, n)
  return(out)
}

getQ = function(K, res, s2){    
  Q = 1/s2*res %*% K %*% res
}

getLambda_davies = function(K, P0){
  PKP = P0 %*% K %*% P0
  ee = eigen(PKP, symmetric = T)         
  lambda0 = ee$values[ee$values > 1e-10]
  return(lambda0)    
}

getindivP_davies = function(Q, lambda0, n, px){
  
  if (length(lambda0) >= n-px){ 
    # In fact, it can not go bigger than n-p because P0 has rank n-p
    lambda = c(lambda0 - Q/(n-px))
    k = length(lambda)
    h = rep(1,k)
  }else{
    lambda = c(lambda0 - Q/(n-px), -Q/(n-px))
    k = length(lambda0)
    h = c(rep(1, k), n-px - k)
  }
  
  p_davies = davies(q = 0, lambda = lambda, h = h, acc = 1e-6)$Qq  
  p_davies = ifelse(p_davies < 0, 0, p_davies) 
  
  return(p_davies)  
}

getSat = function(Q, keppa_tlt, niu_tlt){
  p_sat= 1 - pchisq(Q/keppa_tlt, df = niu_tlt)
  p_sat = ifelse(p_sat<0, 0, p_sat)
  return(p_sat)
}

getParamSat = function(K, P0,px){
  n        = nrow(K)
  POK      = P0%*%K 
  e_tlt    = sum(diag(POK))/2
  Iss      = 0.5*(n-px)#.5*sum(diag(P0 %*% P0))
  W        = POK %*% P0
  Its      = .5*sum(diag(W))
  Itt      = .5*sum(diag(W %*% K))
  Itt_tlt  = Itt - Its^2/Iss
  niu_tlt  = 2*e_tlt^2/Itt_tlt
  keppa_tlt= Itt_tlt/e_tlt
  return(list(niu_tlt = niu_tlt, keppa_tlt = keppa_tlt))
}
# This is to permute the residuals 
getQsim_continuous = function(mod,  nperm, X1, Ks){
  res = resid(mod)
  n = length(res)
  perm    = sapply(1:nperm, function(x) permuted.index(n))
  y_star  = mod$fitted + matrix(res[perm],n,nperm)
  
  res_sim = qr.resid(qr(X1), y_star) 
  px = NCOL(X1)
  modelVar= function(x, px){sum((x - mean(x))^2)/(n - px)}  
  
  sigma2_sim= apply(res_sim, 2, modelVar, px) # Already the variance
  
  Q_sim   = sapply(1:length(Ks), function(j){
    sapply(1:nperm, function(i){
      res_sim[,i] %*% Ks[[j]] %*%res_sim[,i]/sigma2_sim[i]})
  })
  return(Q_sim)    
}

# For binary outcome 
Get_Var_Elements =function(m4,u1,u2){
  temp1 = u1^2 * u2^2
  a1    = sum(m4 * temp1)
  a2    = sum(u1^2) * sum(u2^2) - sum(temp1)
  a3    = sum(u1*u2)^2 - sum(temp1)
  return(a1+a2+2*a3)
}
getHm = function(Q,muQ, varQ, df){
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p = 1 - pchisq(Q_corrected ,df = df)
  p = ifelse(p<0, 0, p)
  return(p)
}
getIndivP_hm = function(K, res, mu, D0, P0){
  Q   = t(res)%*% K %*% res
  # K1  = 1/D0 * P01 %*% K %*% t(1/D0 * P01 )
  K1 = P0 %*% (D0*t(D0*K)) %*% P0
  eK  = eigen(K1, symmetric = T)
  # Instead of matching the first two moments, match to the fourth moment
  # Code adapted from SKAT package
  lambda = eK$values[eK$values > 1e-10]
  U   = as.matrix(eK$vectors[,eK$values > 1e-10])
  p.m = length(lambda)
  m4  = (3*mu^2-3*mu +1)/(mu*(1-mu))
  
  zeta =rep(0,p.m)
  var_i=rep(0,p.m)
  varQ = 0
  
  for(i in 1:p.m){   # The diagonals
    temp.M1 = sum(U[,i]^2)^2 - sum(U[,i]^4)
    zeta[i] = sum(m4 * U[,i]^4) + 3* temp.M1 # because ( \sum .)^4, not ^2
    var_i[i]= zeta[i] - 1
  }
  
  if(p.m == 1){
    Cov_Mat = matrix(zeta* lambda^2, ncol=1,nrow=1)
  } else if(p.m > 1){
    Cov_Mat = diag(zeta* lambda^2)
    for(i in 1:(p.m-1)){
      for(j in (i+1):p.m){
        Cov_Mat[i,j] = Get_Var_Elements(m4,U[,i],U[,j])
        Cov_Mat[i,j] = Cov_Mat[i,j]* lambda[i]* lambda[j]
      }
    }
  }
  Cov_Mat       = Cov_Mat + t(Cov_Mat)
  diag(Cov_Mat) = diag(Cov_Mat)/2
  
  varQ = sum(Cov_Mat) - sum(lambda)^2
  muQ  = sum(lambda)
  lambda.new = lambda * sqrt(var_i)/sqrt(2)
  df         =  sum(lambda.new^2)^2/sum(lambda.new^4)
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p_corrected= 1 - pchisq(Q_corrected ,df = df)
  
  p_corrected = ifelse(p_corrected <0, 0, p_corrected)
  return(list(p_hm= p_corrected, Q = Q, muQ = muQ, varQ = varQ, df = df))
}

getIndivP_binary = function(K, res, D0, px, P0){
  n = length(res)
  Q <- as.numeric(res %*% K %*% res) 
  PKP = P0 %*% (D0*t(D0 * K)) %*% P0 
  eP0 = c(rep(1, n-px), rep(0, px))
  ePKP = eigen(PKP, symmetric = T)$values
  lambda0 = ePKP - Q*eP0/n  # the MLE of s2 is 1, therefore, we should divide by .        
  lambda0 = lambda0[abs(lambda0) >= 1e-10]
  p = davies(0, lambda=lambda0, acc = 1e-6)$Qq
  p = ifelse(p < 0, 0, p)
  return(list(Q = Q, ePKP = ePKP, p = p ))
}
