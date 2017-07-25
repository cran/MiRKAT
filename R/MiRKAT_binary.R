
MiRKAT_binary = function(y, X = NULL, Ks, family= "binomial", nperm = 999, method = "davies"){
  n <- length(y)
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
  
  
  qX1 <- qr(X1)
  ## Take care of aliased variables and pivoting in rhs
  X1 <- X1[, qX1$pivot, drop=FALSE]
  X1 <- X1[, 1:qX1$rank, drop=FALSE]
  options(warn=2)  # make sure this model is correct
  mod <- glm(y ~ X1-1, family = binomial)
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
      return(indivP = ps)
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

    return(list(indivP = ps , omnibus_p = p_final))    
  }
  if (method == "moment"){
    
    S = sapply(Ks, getIndivP_hm, res, mu, D0, P0)
    ps = as.numeric(unlist(S[1,]))
    Qs = unlist(S[2,])
    
    if (length(Ks) == 1){
      return(indivP = ps)
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

    return(list(indivP = ps, omnibus_p = p_final))
  }

  if (method == "permutation"){
    
    Qs = lapply(Ks, getQ,res, s2 = 1)
    
    q_sim = sapply(1:nperm, function(i){
      ind <- sample(n)
      p1 = sapply(1:length(Ks), function(j) {
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
        return(Q1)
      })  
    }) 

    if (length(Ks) == 1){
     p_perm = (sum(q_sim > Qs)+ 1)/(nperm + 1) 
     return(indivP = p_perm)
    }
  
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)   
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm + 1)
    
    return(list(indivP = p_perm , omnibus_p = p_final)) 
    
  }
}

