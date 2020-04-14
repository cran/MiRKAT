#'D2K
#'
#'Construct kernel matrix from distance matrix.
#'
#'Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. The kernel matrix is constructed as \eqn{K = -(I-11'/n)D^2(I-11'/n)/2}, where D is the pairwise distance matrix, I is the identity
#'matrix, and 1 is a vector of ones.
#'
#'\eqn{D^2} represents element-wise square.
#'
#'To ensure that \eqn{K} is positive semi-definite, a positive semi-definiteness correction is conducted
#'
#'@param D An n by n matrix giving pairwise distances or dissimilarites, where n is sample size.
#'
#'@return An n by n kernel or similarity matrix corresponding to the distance matrix given.
#'
#'@author 
#'Ni Zhao
#'
#'@references
#'Zhao, Ni, et al. "Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test
#'
#'@examples
#'library(GUniFrac)
#'
#'#Load in data and create a distance matrix
#'data(throat.tree)
#'data(throat.otu.tab)
#'unifracs <- GUniFrac(throat.otu.tab, throat.tree, alpha=c(1))$unifracs
#'D1 <- unifracs[,,"d_1"]
#'
#'#Function call
#'K <- D2K(D1)
#'
#'@export
#'
D2K <-
function(D){
  n <- nrow(D)
  centerM <- diag(n) - 1/n
  K <- -0.5*centerM %*% (D*D) %*% centerM
  eK <- eigen(K, symmetric=TRUE)
  K <- eK$vector %*% diag(pmax(0,eK$values)) %*% t(eK$vector)
  return(K)
}
