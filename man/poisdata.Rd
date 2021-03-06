% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poisdata.R
\docType{data}
\name{poisdata}
\alias{poisdata}
\title{Simulated DEPENDENT data with POISSON (count) traits for correlated regression-based analysis (i.e. CSKAT, GLMMMiRKAT)}
\format{
A list containing three data objects for correlated microbiome data with binary response variable (described below). 

\describe{
  \item{pois.otu.tab}{Simulated OTU data for correlated regression-based analysis; 59 rows and 730 columns, rows being patients and 
  columns being OTUs}
  \item{pois.meta}{Simulated metadata for correlated regression-based analysis; 59 rows and 4 columns, rows being patients and 
  columns being the outcome variable, subject identifier, and covariates to possibly account for in any regression modeling}
  \item{pois.tree}{Simulated rooted phylogenetic tree with 730 tips and 729 nodes}
}
}
\usage{
data(poisdata)
}
\description{
Simulated DEPENDENT data with POISSON (count) traits for correlated regression-based analysis (i.e. CSKAT, GLMMMiRKAT)
}
\keyword{datasets}
