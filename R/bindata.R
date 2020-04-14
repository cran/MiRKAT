#' Simulated DEPENDENT data with BINOMIAL traits for correlated regression-based analysis (i.e. CSKAT, GLMMMiRKAT)
#'
#' @usage data(bindata)
#' 
#' @format A list containing three data objects for correlated microbiome data with binary response variable (described below). 
#' 
#'\describe{
#'   \item{bin.otu.tab}{Simulated OTU data for correlated regression-based analysis; 59 rows and 730 columns, rows being patients and 
#'   columns being OTUs}
#'   \item{bin.meta}{Simulated metadata for correlated regression-based analysis; 59 rows and 4 columns, rows being patients and 
#'   columns being the outcome variable, subject identifier, and covariates to possibly account for in any regression modeling}
#'   \item{bin.tree}{Simulated rooted phylogenetic tree with 730 tips and 729 nodes}
#'}
#'
#'
#'
"bindata"