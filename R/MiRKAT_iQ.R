#'@title MiRKAT-iQ
#'@description Integrated quantile regression-based kernel association test.
#'
#'@param Y A numeric vector of the continuous outcome variable.
#'@param X A numeric matrix for additional covariates that you want to adjust for.
#'@param K A list of n by n kernel matrices at a single n by n kernel matrix, where n is the sample size.
#'@param weight A length 4 vector specifying the weight for Cauchy combination, corresponding to wilcoxon/normal/inverselehmann/lehmann functions. The sum of the weight should be 1.
#'
#'@return Returns a list containing the p values for single kernels, or the omnibus p-value if multiple candidate kernel matrices are provided.
#'@author 
#'Tianying Wang, Xiang Zhan. 
#' 
#' @references
#'Wang T, et al. (2021) Testing microbiome association using integrated quantile regression models. Bioinformatics (to appear). 
#' 
#'@import quantreg
#'@examples
#'
#'library(GUniFrac)
#'library(quantreg)
#'library(PearsonDS)
#'library(MiRKAT)
#'
#'data(throat.tree)
#'data(throat.otu.tab)
#'
#'## Create UniFrac and Bray-Curtis distance matrices 
#'unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
#'if (requireNamespace("vegan")) {
#'  library(vegan)
#'  BC= as.matrix(vegdist(throat.otu.tab, method="bray"))
#'  Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], BC = BC)
#'} else {
#'  Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"])
#'}
#'
#'## Convert to kernels           
#'Ks = lapply(Ds, FUN = function(d) D2K(d))
#'covar = cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male"))
#'n = nrow(throat.meta)
#'y = rnorm(n)
#'result = MiRKAT.iQ(y, X = covar, K = Ks)
#'
#'@export
#'
MiRKAT.iQ = function(Y, X, K, weight = c(0.25, 0.25, 0.25, 0.25)){
  if(is.list(K)){
    nk = length(K)
    p.val.singleK = matrix(0, nrow = nk, ncol = 5)
    for(ik in 1:nk){
      results.temp = rep(0,5)
      results.temp[1] = fKMQR.test(Y,X,tau = -1, score = "wilcoxon", K = K[[ik]])
      results.temp[2] = fKMQR.test(Y,X,tau = -1, score = "normal", K = K[[ik]])
      results.temp[3] = fKMQR.test(Y,X,tau = -1, score = "inverselehmann", K = K[[ik]])
      results.temp[4] = fKMQR.test(Y,X,tau = -1, score = "lehmann", K = K[[ik]])
      
      results.temp[5] = pcauchy(sum(weight*(tan((0.5-results.temp[1:4])*pi))), lower.tail = FALSE)
      
      
      p.val.singleK[ik, ] = results.temp
      
    }
    colnames(p.val.singleK) = c("wilcoxon", "normal", "inverselehmann", "lehmann", "iQ")
    rownames(p.val.singleK) = names(K)
    p.val = as.numeric(apply(p.val.singleK, 2, function(x){pcauchy(mean(tan((0.5-x)*pi)), lower.tail = FALSE)}))
    names(p.val) = c("wilcoxon", "normal", "inverselehmann", "lehmann", "iQ")
    p.list = list(p.val = p.val, p.val.singleK = p.val.singleK)
  }else{
    results.temp = rep(0,5)
    results.temp[1] = fKMQR.test(Y,X,tau = -1, score = "wilcoxon", K = K)
    results.temp[2] = fKMQR.test(Y,X,tau = -1, score = "normal", K = K)
    results.temp[3] = fKMQR.test(Y,X,tau = -1, score = "inverselehmann", K = K)
    results.temp[4] = fKMQR.test(Y,X,tau = -1, score = "lehmann", K = K)
    
    results.temp[5] = pcauchy(sum(weight*(tan((0.5-results.temp[1:4])*pi))), lower.tail = FALSE)
    
    
    names(results.temp) = c("wilcoxon", "normal", "inverselehmann", "lehmann", "iQ")
    
    p.list = list(p.val = results.temp)
  }
  
  return(p.list)
}


fKMQR.test = function(Y, Z, tau, score = NULL, K = NULL) {
  ## fast Kernel machine quantile regression
  ### X nonparamteric var (nxp)
  ## Z parametric var (nxq)
  ## Y response var (nx1)
  ## tau quantile
  
  ## Define some auxiliary functions: tr and KRV function
  KRV=function(K,L){
    n=nrow(K)
    A=scale(K,scale=F) ## that is A=HK
    W=scale(L,scale=F) ## that is W=HL
    Fstar=tr(A%*%W)
    mean.krv=tr(A)*tr(W)/(n-1)	## mean of KRV
    T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
    Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
    temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
    temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
    temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
    temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
    temp2=temp21*temp22/temp23
    variance.krv=temp1+temp2		## variance of KRV
    T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
    T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
    t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
    t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
    t3=24*(n^2-n-4)*(U*Bs+B*Us)
    t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
    t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
    t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
    t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
    t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
    t8=24*(t81+t82)
    t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
    t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
    t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
    t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
    t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
    t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
    t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
    t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
    t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
    t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
    t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
    t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
    t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
    t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
    t20=-(n-2)*(t201+t202+t203)
    temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
    temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
    mom3=temp31/temp32
    skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
    m1=mean.krv
    m2=variance.krv
    m3=skewness.krv
    shape=4/m3^2
    scale=sqrt(m2)*m3/2
    location=m1-2*sqrt(m2)/m3
    PIIIpars=list(shape,location,scale)
    pv=1-ppearsonIII(Fstar, params=PIIIpars)
    return(pv)
  }
  
  tr=function(x){return(sum(diag(x))) }
  
  n=length(Y)
  #if(!is.null(K)){
  #  K = IBS.kernel(X) ##kernel matrix
  #}
  
  if(tau == -1){
    null.fit = rq(Y~Z, tau = -1)
    null.fit = transfer_v(null.fit)
    w = NewRanks(null.fit, score = score)
    w = w$ranks/sqrt(w$A2)
  }else{
    quantresult<-rq(Y~Z, tau=tau)
    parametersolve<-quantresult$coefficients
    residual<-quantresult$residuals
    
    w<-rep(0,n)
    for (i in 1:n) {
      if (residual[i]==0)
        residual[i]<-rbinom(1,1,1-tau)-0.5
      w[i]<-tau-(residual[i]<0)
    }
  }
  
  Kw=w%*%t(w)
  pv=KRV(Kw,K)
  return(c(pv))
}




