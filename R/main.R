
MiRKAT = function(y, X = NULL, Ks, out_type = "C", nperm = 999, method = "davies"){
  n = length(y)
  
  if (any(is.na(y))){
    ids = which(is.na(y))
    stop(paste("subjects", ids, "has missing response, please remove before proceed \n")) 
  }
  
  
  if(is.null(X)==FALSE){
    if(NROW(X)!= length(y)) stop("Dimensions of X and y don't match.")
  }
  
  if (class(Ks) == "matrix"){
    Ks = list(Ks)
  }
  
  if (class(Ks) == "list"){  
    if((any(lapply(Ks, "nrow")!= n))|(any(lapply(Ks,  "ncol")!= n))){
      stop("distance matrix need to be n x n, where n is the sample size \n ")    
    } 
    if (class(Ks) != "list"){
      stop("Distance needs to be a list of n x n matrices or a single n x n matrix \n")  
    }
  }
  

  if (!is.null(X)){
    if (any(is.na(X))){
      stop("NAs in  covariates X, please impute or remove subjects which has missing covariates values") 
    }  
  }
  
  if (method == "moment" & n < 100 & out_type == "C"){
    
    warning("Continuous outcome: sample size < 100, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  if (method == "moment" & n < 200 & out_type == "D"){
    
    warning("Continuous outcome: sample size < 200, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  
  if (!(out_type %in%  c("C", "D"))){
    stop("Currently only continuous and Binary outcome are supported. Please choose out_type = \"C\" or \"D\" ")
  }
  if(out_type  == "C"){
    re = MiRKAT_continuous(y, X = X, Ks = Ks, method = method, nperm = nperm)  
  }
  
  if(out_type  == "D"){
    re = MiRKAT_binary(y, X = X, Ks = Ks, method = method, nperm = nperm)  
  }
  
  return(re)
}
