\name{MMiRKAT}
\alias{MMiRKAT}

\title{Multivariate Microbiome Regression-based Kernel Association Test}
\description{
Test for association between overall microbiome composition and multiple continuous outcomes.
}
\usage{
MMiRKAT(Y, X=NULL, K)
}

\arguments{
    \item{Y}{A numerical n by p matrix of p continuous outcome variables.}
    \item{X}{A numerical n by q matrix or data frame, containing q additional covariates that you want to adjust for (Default = NULL). If it is NULL, a intercept only model was fit.}
    \item{K}{A numerical n by n kernel matrix, where n is the sample size. It can be constructed from microbiome data base on Bray-Curtis or UniFrac distance of microbiome data.}
}
\details{
   Y and X (if not NULL) should all be numerical matrices with the same number of rows n.

   K should be a single n by n matrix. If you have a distance matrix from metagenomic data, a kernel can be constructed through function D2K. The kernel can also be constructed through other mathematical approaches.

   Missing data is not permitted. Please remove all individuals with missing Y, X, K prior to analysis

   The method of generating kernel specific p-values is "davies", which represents an exact method that computes the p-value by inverting the characteristic function of the mixture chisq.
}
\value{
P-value obtained using small-sample correction
}
\references{
  Zhan, X., Tong, X., Zhao, N., Maity,A., Wu, M.C., and Chen,J. A small-sample multivariate kernel machine test for microbiome association studies. Genetic Epidemiology, 41(3), 210-220. DOI: 10.1002/gepi.22030
}
\author{
Haotian Zheng, Xiang Zhan, Ni Zhao
}
\examples{
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

Y = matrix(rnorm(n * 3, 0, 1),n ,3)

#############################################################

MMiRKAT(Y = Y, K = K.weighted, X = cbind(Male, anti)) 
}

