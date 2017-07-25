\name{KRV}
\alias{KRV}

\title{
Kernel RV Coefficient Test
}
\description{
kernel RV coefficient test to evaluate the overall association between microbiome composition and high-dimensional or structured phenotype.
}
\usage{
KRV(kernel.otu, y = NULL, X = NULL, kernel.y)
}

\arguments{
  \item{kernel.otu}{
  A numerical n by n kernel matrix. It can be constructed from microbiome data, such as by transforming from a distance metric.
}
  \item{y}{
  A numerical n by p matrix of p continuous phenotype variables (Default = NULL). If it is NULL, a phenotype kernel matrix must be entered for "kernel.y". No need to provide if kernel.y is a matrix.
}
  \item{X}{
  A numerical n by q matrix, containing q additional covariates that you want to adjust for (Default = NULL). If it is NULL, a intercept only model was fit. Covariates can't be adjusted for if kernel.y is a matrix.
  }
  \item{kernel.y}{
  Either a numerical n by n kernel matrix of phenotype or a method to compute the kernel of phenotype.
  Gaussian kernel (kernel.y="Gaussian") can capture general relationship between microbiome and phenotypes; and linear kernel (kernel.y="linear") can be preferred if the underlying relationship is close to linear.
}
}
\details{
   kernel.otu should be a numerical n by n kernel matrix, where n is sample size.  

   When kernel.y is a method ("Gaussian" or "linear") to compute the kernel of phenotype, y should be a numerical phenotype matrix, and X (if not NULL) should be a numerical matrix of covariates. Both y and X should have n rows. 
   
   When kernel.y is a kernel matrix of phenotype, there is no need to provide X and y, and they will be ignored if provided. In this case, kernel.y and kernel.otu should both be numerical matrices with the same number of rows and columns. 

   Missing data is not permitted. Please remove all individuals with missing kernel.otu, y (if not NULL), X (if not NULL), and kernel.y (if a matrix is entered) prior to analysis.
}
\value{
P-value calculated from approximated Pearson type III density
}
\references{
  Zhan, X., Plantinga, A., Zhao, N., and Wu, M.C. A Fast Small-Sample Kernel Independence Test for Microbiome Community-Level Association Analysis. Biometrics. 2017 Mar 10. doi: 10.1111/biom.12684.
}
\author{
Haotian Zheng, Xiang Zhan, Ni Zhao
}
\examples{
library(MASS)
library(GUniFrac)
data(throat.tree)
data(throat.otu.tab)
data(throat.meta)
attach(throat.meta)

set.seed(123)
n = nrow(throat.otu.tab)
Male = (Sex == "Male")**2
Smoker =(SmokingStatus == "Smoker") **2
anti =  (AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")^2
cova = cbind(Male,  anti)

otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff
unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs

D.weighted = unifracs[,,"d_1"]
D.unweighted = unifracs[,,"d_UW"]
D.BC= as.matrix(vegdist(otu.tab.rff , method="bray"))

K.weighted = D2K(D.weighted)
K.unweighted = D2K(D.unweighted)
K.BC = D2K(D.BC)

rho = 0.2
Va = matrix(rep(rho, (2*n)^2), 2*n, 2*n)+diag(1-rho, 2*n)
G = mvrnorm(n, rep(0, 2*n), Va)

#############################################################

KRV(kernel.otu = K.weighted, y = G, X = cova, kernel.y = "Gaussian") 
KRV(kernel.otu = K.weighted, kernel.y = G \%*\% t(G)) 
}