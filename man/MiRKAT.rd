 \name{MiRKAT}
 \alias{MiRKAT}
 \title{Microbiome Regression-based Kernel Association Test}
 \description{
     Test for association between microbiome composition and a continuous or dichotomous outcome by incorporating phylogenetic or nonphylogenetic distance between different microbiomes.
 }
 \references{

  Zhao, N., Chen, J.,Carroll, I. M., Ringel-Kulka, T., Epstein, M.P., Zhou, H., Zhou, J. J.,  Ringel, Y., Li, H. and Wu, M.C. (2015)). Microbiome Regression-based Kernel Association Test (MiRKAT). American Journal of Human Genetics, 96(5):797-807

  Chen, J., Chen, W., Zhao, N., Wu, M~C.and Schaid, D~J. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome Association Studies. 40: 5-19. doi: 10.1002/gepi.21934

Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, \emph{ Journal of the Royal Statistical Society. Series C }, 29, 323-333.

Satterthwaite, F. (1946). An approximate distribution of estimates of variance components. Biom. Bull. 2, 110-114.

Lee S,  Emond MJ,  Bamshad MJ, Barnes KC, Rieder MJ, Nickerson DA; NHLBI GO Exome Sequencing Project-ESP Lung Project Team, Christiani DC, Wurfel MM, Lin X. (2012) Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. \emph{American Journal of Human Genetics}, 91, 224-237.

Zhou, J. J. and Zhou, H.(2015) Powerful Exact Variance Component Tests for the Small Sample Next Generation Sequencing Studies (eVCTest), in submission.

}
\usage{
 MiRKAT(y, X = NULL, Ks, out_type = "C", nperm = 999, method = "davies")
}
\arguments{
    \item{y}{A numeric vector of the a continuous or dichotomous outcome variable.}
    \item{X}{A numerical matrix or data frame, containing additional covariates that you want to adjust for (Default = NULL). If it is NULL, a intercept only model was fit.}
    \item{Ks}{a list of n by n kernel matrices (or a single n by n kernel matrix), where n is the sample size. It can be constructed from microbiome data through distance metric or other approaches, such as linear kernels or Gaussian kernels.}
    \item{out_type}{an indicator of the outcome type. "C" for the continuous outcome and "D" for the dichotomous outcome.}
    \item{nperm}{the number of permutations if method = "permutation" or when multiple kernels are considered. if method = "davies" or "moment", nperm is ignored.}
    \item{method}{a method to compute the kernel specific p-value (Default= "davies").
      "davies" represents an exact method that  computes the p-value by inverting the characteristic function of the mixture chisq. We adopt an exact variance component tests because most of the studies concerning microbiome compositions have modest sample size.
      "moment" represents an approximation method that matches the first two moments.
      "permutation" represents a permutation approach for p-value calculation.}
}

\details{
   y and X (if not NULL) should all be numerical matrices or vectors with the same number of rows.

   Ks should be a list of n by n matrices or a single matrix. If you have distance metric from metagenomic data, each kernel can be constructed through function D2K.  Each kernel can also be constructed through other mathematical approaches.

   Missing data is not permitted. Please remove all individuals with missing y, X, Ks prior to analysis

   Parameter "method" only concerns with how kernel specific p-values are generated. When Ks is a list of multiple kernels, omnibus p-value is computed through permutation from each individual p-values, which are calculated through method of choice.
}

\value{
\item{indivP}{p-value from each candidate kernel}
\item{omnibus_p}{omnibus p value by considering multiple candidate kernels.}
}

\author{Ni Zhao}

\examples{

#############################################################
# Generate data
set.seed(1)
n = 100
family= "binomial"
nperm = 999
X = rnorm(n)
Z = matrix((runif(n*5) > 0.5)^2, n, 5)
K = Z \%*\% t(Z)
K2 = (Z \%*\% t(Z)+ 1)^2
Ks = list(K, K2)
y = rnorm(n)

MiRKAT(y, X = X, Ks = Ks, out_type="C", method = "davies")

#################################
y = (runif(n) < 0.5)^2
MiRKAT(y, X = X, Ks = Ks, out_type="D")
}


