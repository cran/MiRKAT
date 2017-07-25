\name{D2K}
\alias{D2K}
\title{Construct kernel matrix from distance metric}
\description{Construct kernel matrix from distance matrix (matrix of pairwise distances) for microbiome data.}
\usage{D2K(D)}
\arguments{
  \item{D}{An n by n matrix giving pairwise distances or dissimilarities, where n is the sample size.}
}
\details{
The kernel matrix is constructed as \eqn{K = - (I - 11'/n)D^2 (I - 11'/n)/2}, where D is the pairwise distance matrix, I is an identity matrix and 1 is a vector of 1. \eqn{D^2} represents element wise square. To ensure that \eqn{K} to be positive semi-definite, a positive semi-definiteness correction is conducted that the negative eigen values of \eqn{K} are replaced by zeros.

}
\value{An n by n kernel or similarity matrix corresponding to the distance matrix.}
\references{
Zhao, Ni, et al. "Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test." The American Journal of Human Genetics 96.5 (2015): 797-807.
}
\author{Ni Zhao}

\examples{
  #############################################################
 require(GUniFrac)
 # Load in data and create a distance metric
 data(throat.tree)
 data(throat.otu.tab)
 unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
 D1 = unifracs[,,"d_1"]
 # Function
 K = D2K(D1)
}