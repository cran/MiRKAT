% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiRKAT.R
\name{MiRKAT}
\alias{MiRKAT}
\title{Microbiome Regression-based Kernel Association Test}
\usage{
MiRKAT(
  y,
  X = NULL,
  Ks,
  out_type = "C",
  method = "davies",
  omnibus = "permutation",
  nperm = 999,
  returnKRV = FALSE,
  returnR2 = FALSE
)
}
\arguments{
\item{y}{A numeric vector of the a continuous or dichotomous outcome variable.}

\item{X}{A numeric matrix or data frame, containing additional covariates that you want to adjust for. If NULL, a intercept
only model is used. Defaults to NULL.}

\item{Ks}{A list of n by n kernel matrices or a single n by n kernel matrix, where n is the sample size. It can be constructed 
from microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.}

\item{out_type}{An indicator of the outcome type ("C" for continuous, "D" for dichotomous).}

\item{method}{Method used to compute the kernel specific p-value. "davies" represents an exact method that computes the p-value by
inverting the characteristic function of the mixture chisq. We adopt an exact variance component tests because most of the studies
concerning microbiome compositions have modest sample size. "moment" represents an approximation method that matches the first
two moments. "permutation" represents a permutation approach for p-value calculation. Defaults to "davies".}

\item{omnibus}{A string equal to either "cauchy" or "permutation" (or nonambiguous abbreviations thereof), specifying whether 
to use the Cauchy combination test or residual permutation to generate the omnibus p-value.}

\item{nperm}{The number of permutations if method = "permutation" or when multiple kernels are considered. If method = "davies" or 
"moment", nperm is ignored. Defaults to 999.}

\item{returnKRV}{A logical indicating whether to return the KRV statistic (a measure of effect size). Defaults to FALSE.}

\item{returnR2}{A logical indicating whether to return R-squared. Defaults to FALSE.}
}
\value{
Returns a list containing the following elements: 
   \item{p_values}{P-value for each candidate kernel matrix}
   \item{omnibus_p}{Omnibus p-value considering multiple candidate kernel matrices}
   \item{KRV}{Kernel RV statistic (a measure of effect size). Only returned if returnKRV = TRUE.}
   \item{R2}{R-squared. Only returned if returnR2 = TRUE.}
}
\description{
Test for association between microbiome composition and a continuous or dichotomous outcome by incorporating 
phylogenetic or nonphylogenetic distance between different microbiomes.
}
\details{
y and X (if not NULL) should all be numeric matrices or vectors with the same number of rows.

Ks should be a list of n by n matrices or a single matrix. If you have distance metric(s) from metagenomic data, each kernel can be
constructed through function D2K. Each kernel can also be constructed through other mathematical approaches.

Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis

Parameter "method" only concerns with how kernel specific p-values are generated. When Ks is a list of multiple kernels, omnibus 
p-value is computed through permutation from each individual p-values, which are calculated through method of choice.
}
\examples{
library(GUniFrac)

data(throat.tree)
data(throat.otu.tab)
data(throat.meta)

unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
if (requireNamespace("vegan")) {
 library(vegan)
 BC= as.matrix(vegdist(throat.otu.tab, method="bray"))
 Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], BC = BC)
} else {
 Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"])
}
Ks = lapply(Ds, FUN = function(d) D2K(d))

covar = cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male"))

# Continuous phenotype
n = nrow(throat.meta)
y = rnorm(n)
MiRKAT(y, X = covar, Ks = Ks, out_type="C", method = "davies")

# Binary phenotype 
y = as.numeric(runif(n) < 0.5)
MiRKAT(y, X = covar, Ks = Ks, out_type="D")


}
\references{
Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J., Ringel, Y., Li, H. and Wu, M.C. (2015)).
Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807

Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome
Association Studies. 40: 5-19. doi: 10.1002/gepi.21934

Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal 
Statistical Society. Series C , 29, 323-333.

Satterthwaite, F. (1946). An approximate distribution of estimates of variance components. Biom. Bull. 2, 110-114.

Lee S, Emond MJ, Bamshad MJ, Barnes KC, Rieder MJ, Nickerson DA; NHLBI GO Exome Sequencing Project-ESP Lung Project Team, 
Christiani DC, Wurfel MM, Lin X. (2012) Optimal unified approach for rare variant association testing with application to small 
sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.

Zhou, J. J. and Zhou, H.(2015) Powerful Exact Variance Component Tests for the Small Sample Next Generation Sequencing Studies 
(eVCTest), in submission.
}
\author{
Ni Zhao
}
