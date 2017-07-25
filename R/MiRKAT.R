

MiRKAT_continuous = function(y, X = NULL, Ks,method = "davies", nperm = 999){
  n = length(y) 

  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
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
  
  if (class(Ks) == "matrix"){Ks = list(Ks)}
  Qs = lapply(Ks, getQ, res, s2)  
  if (method == "davies"){    
    lambda0 = lapply(Ks, getLambda_davies, P0)
    if (length(Ks) == 1){
      p = getindivP_davies(Qs[[1]], lambda0[[1]], n, px)    
      return(indivP = p)  # This should prevent the algorithm going any further
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
   
    return(list(indivP = ps , omnibus_p = p_final))    
  } 
  if (method == "moment"){ 
    parm    =  sapply(Ks, getParamSat, P0, px)  
    if (length(Ks) == 1){
      p_sat = getSat(Qs[[1]], keppa_tlt = parm[2,1]$keppa_tlt, niu_tlt = parm[1,1]$niu_tlt)  
      return(indivP =as.numeric(unlist(p_sat)))
     }
    p_sat = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      p_sat[i] = getSat(Qs[[i]], keppa_tlt = parm[2,i]$keppa_tlt, niu_tlt = parm[1,i]$niu_tlt)    
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
    
    return(list(indivP = p_sat , omnibus_p = p_final)) 
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
      return(indivP = p_perm)    
      
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
  
  return(list(indivP = p_perm , omnibus_p = p_final)) 

  }
  
}
