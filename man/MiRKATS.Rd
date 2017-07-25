\name{MiRKATS}
\alias{MiRKATS}
\title{Microbiome Regression-based Kernel Association Test for Survival}
\description{
Community level test for association between microbiome composition and survival outcomes (right-censored time-to-event data) using kernel matrices to compare similarity between microbiome profiles with similarity in survival times.
}
\usage{
MiRKATS(kd, distance = FALSE, obstime, delta, covar = NULL, beta = NULL,
perm = FALSE, nperm = 1000)
}
\arguments{
  \item{kd}{A numeric n by n kernel matrix or matrix of pairwise distances/dissimilarities (where n is the sample size).}
  \item{distance}{Logical, indicating whether kd is a distance matrix (default = FALSE).}
  \item{obstime}{A numeric vector of follow-up (survival/censoring) times.}
  \item{delta}{Event indicator: a vector of 0/1, where 1 indicates that the event was observed for a subject (so "obstime" is survival time), and 0 indicates that the subject was censored.}
  \item{covar}{A vector or matrix of numeric covariates, if applicable (default = NULL).}
  \item{beta}{A vector of coefficients associated with covariates. If beta is NULL and covariates are present, coxph is used to calculate coefficients (default = NULL).}
  \item{perm}{Logical, indicating whether permutation should be used instead of analytic p-value calculation (default=FALSE). Not recommended for sample sizes of 100 or more.}
  \item{nperm}{Integer, number of permutations used to calculate p-value if perm==TRUE (default=1000).}
}
\details{
obstime, delta, and covar should all have n rows, and the kernel or distance matrix should be a single n by n matrix. If a distance matrix is entered (so distance=TRUE), a kernel matrix will be constructed from the distance matrix.

Missing data is not permitted. Please remove individuals with missing data on y, X or in the kernel or distance matrix prior to using the function.

The Efron approximation is used for tied survival times.
}
\value{P-value obtained using small sample correction}

\references{
  Plantinga, A., Zhan, X., Zhao, N., Chen, J., Jenq, R., and Wu, M.C. MiRKAT-S: a distance-based test of association between microbiome composition and survival times. Microbiome, 2017:5-17. doi: 10.1186/s40168-017-0239-9

  Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J.,  Ringel, Y., Li, H. and Wu, M.C. (2015)). Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807

  Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome Association Studies. 40:5-19. doi: 10.1002/gepi.21934

  Efron, B. (1977) "The efficiency of Cox's likelihood function for censored data." Journal of the American statistical Association 72(359):557-565.

  Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal Statistical Society Series C, 29:323-333.
}

\author{Anna Plantinga}

\examples{

###################################
# Generate data
require(GUniFrac)
set.seed(1)

# Throat microbiome data
data(throat.tree)
data(throat.otu.tab)
unifracs = GUniFrac(throat.otu.tab, throat.tree, alpha = c(1))$unifracs
D1 = unifracs[,,"d_1"]  # 60 subjects

# Covariates and outcomes
X <- matrix(rnorm(120), nrow=60)
S <- rexp(60, 3)
C <- rexp(60, 1)
D <- (S<=C)        # event indicator
U <- pmin(S, C)    # observed follow-up time

MiRKATS(kd = D1, distance = TRUE, obstime = U, delta = D, covar = X, beta = NULL)

}


