% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSKAT.R
\name{CSKAT}
\alias{CSKAT}
\title{Small-sample SKAT for correlated (continuous) data ('c' stands for 'correlated'). Called within GLMM-MiRKAT.}
\usage{
CSKAT(lmer.obj, Ks)
}
\arguments{
\item{lmer.obj}{A fitted lme4 object (model under H0)}

\item{Ks}{A kernel matrix or list of kernels, quantifying the similarities between samples.}
}
\value{
\describe{
  \item{p.value}{Association p-values}
  \item{Q.adj}{Adjusted score statistics}
}
}
\description{
Compute the adjusted score statistic and p-value
}
\references{
Zhan X, et al. (2018) A small-sample kernel association test for correlated data with application to microbiome association studies.  
Genet Epidemiol.
}
\author{
Nehemiah Wilson, Anna Plantinga, Xiang Zhan, Jun Chen.
}
