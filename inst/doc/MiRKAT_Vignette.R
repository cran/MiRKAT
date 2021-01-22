## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ----load-packages, echo=FALSE------------------------------------------------
library(MiRKAT, quietly = TRUE)
library(GUniFrac, quietly = TRUE)
library(vegan, quietly = TRUE)
library(cluster); library(CompQuadForm); library(dirmult); library(lme4)
library(MASS); library(Matrix); library(survival); library(permute)
library(propr); library(kableExtra)

## -----------------------------------------------------------------------------
data("throat.otu.tab")
data("throat.meta")
data("throat.tree")

## -----------------------------------------------------------------------------
Male <- as.numeric(throat.meta$Sex == "Male")
Smoker <- as.numeric(throat.meta$SmokingStatus == "Smoker")
anti <- as.numeric(throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")
covar <- cbind(Male, anti)

## -----------------------------------------------------------------------------
# Note: may want to rarefy the data to the minimum sequencing depth for unweighted (presence/absence) distances 
## otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff
unifracs <- GUniFrac(throat.otu.tab, throat.tree, alpha = c(0, 0.5, 1))$unifracs
D.weighted = unifracs[,,"d_1"]
D.unweighted = unifracs[,,"d_UW"]
D.generalized = unifracs[,,"d_0.5"]
D.BC = as.matrix(vegdist(throat.otu.tab, method="bray"))

## -----------------------------------------------------------------------------
K.weighted = D2K(D.weighted)
K.unweighted = D2K(D.unweighted)
K.generalized = D2K(D.generalized)
K.BC = D2K(D.BC)

## -----------------------------------------------------------------------------
MiRKAT(y = Smoker, X = covar, Ks = K.weighted, out_type = "D", 
       method = "davies", returnKRV = TRUE, returnR2 = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
mytab = matrix(NA, nrow = 6, ncol = 5) 
mytab[1,] = c("Distance", "Phylogeny?", "Abundance?", "Other notes", "References") 
mytab[2,] = c("Unweighted UniFrac", "Yes", "No" , "", "1") 
mytab[3,] = c("Weighted UniFrac", "Yes", "Yes", "", "2") 
mytab[4,] = c("Generalized UniFrac", "Yes", "(Yes)", "Parameter alpha defines extent to which abundance is taken into account", "3") 
mytab[5,] = c("Jaccard", "No", "No", "1 - (taxa in both)/(taxa in either); typically presence/absence, but can be extended to an abundance-weighted version", "4,5") 
mytab[6,] = c("Bray-Curtis", "No", "Yes", "Similar to Jaccard, but uses counts", "6")

kable(mytab, booktabs=TRUE) %>% 
  column_spec(4, width="20em") 

## -----------------------------------------------------------------------------
Ks = list(K.weighted = K.weighted, K.unweighted = K.unweighted, K.BC = K.BC)
MiRKAT(y = Smoker, Ks = Ks, X = covar, out_type = "D", 
       method = "davies", omnibus = "permutation", 
       returnKRV = TRUE, returnR2 = TRUE)

MiRKAT(y = Smoker, Ks = Ks, X = covar, out_type = "D", 
       method = "davies", omnibus = "cauchy", 
       returnKRV = FALSE, returnR2 = FALSE)

## -----------------------------------------------------------------------------
y <- rnorm(nrow(K.weighted))
MiRKAT(y = y, Ks = Ks, X = covar, out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy", returnKRV = TRUE, returnR2 = TRUE)

## -----------------------------------------------------------------------------
# Simulate outcomes
# Here, outcome is associated with covariates but unassociated with microbiota
# Approximately 33% censoring
SurvTime <- rexp(60, (1 + Smoker + Male))
CensTime <- rexp(60, 0.75)
Delta <- as.numeric(SurvTime <= CensTime )
ObsTime <- pmin(SurvTime, CensTime)

## -----------------------------------------------------------------------------
# Davies' exact method 
MiRKATS(obstime = ObsTime, delta = Delta, X = cbind(Smoker, Male, anti), Ks = Ks, 
        perm = F, omnibus = "cauchy", returnKRV = T, returnR2 = T)

## -----------------------------------------------------------------------------
# Permutation 
MiRKATS(obstime = ObsTime, delta = Delta, X = cbind(Smoker, Male, anti), Ks = Ks, 
        perm = T, omnibus = "cauchy", returnKRV = T, returnR2 = T)

## -----------------------------------------------------------------------------
set.seed(1)
n <- nrow(throat.otu.tab)
Y <- matrix(rnorm(n*3, 0, 1), n, 3)
MMiRKAT(Y = Y, X = cbind(Male, anti), Ks = Ks, 
        returnKRV = TRUE, returnR2 = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)
rho = 0.2
Va = matrix(rep(rho, (2*n)^2), 2*n, 2*n)+diag(1-rho, 2*n)
G = mvrnorm(n, rep(0, 2*n), Va)

## -----------------------------------------------------------------------------
lin.K.y = G %*% t(G)
KRV(kernels.otu = Ks, kernel.y = lin.K.y, 
    omnibus = "kernel_om", 
    returnKRV = TRUE, returnR2 = TRUE)

## -----------------------------------------------------------------------------
KRV(y = G, X = cbind(Male, anti), kernels.otu = Ks, kernel.y = "linear")

## -----------------------------------------------------------------------------
y <- rchisq(nrow(K.weighted), 1)
MiRKAT.R(y = y, X = cbind(Male, anti), Ks = Ks)

## -----------------------------------------------------------------------------
MiRKAT.Q(y = y, X= cbind(Male, anti), Ks = Ks)

## -----------------------------------------------------------------------------
# Import example microbiome data with Gaussian traits
data(nordata)

# Extract microbiome and meta information
otu.tab <- nordata$nor.otu.tab
tree <- nordata$nor.tree
meta <- nordata$nor.meta

## -----------------------------------------------------------------------------
Ds <- GUniFrac(otu.tab, tree, alpha=c(0.5, 1))$unifracs
 
D.unweighted <- Ds[,,"d_UW"]
D.generalized <- Ds[,,"d_0.5"]
D.weighted <- Ds[,,"d_1"]

## -----------------------------------------------------------------------------
K.unweighted <- D2K(D.unweighted)
K.generalized <- D2K(D.generalized)
K.weighted <- D2K(D.weighted)
Ks <- list(K.weighted = K.weighted, 
           K.unweighted = K.unweighted, 
           K.generalized = K.generalized)

## -----------------------------------------------------------------------------
# GLMM-MiRKAT with computational p-values (Davies + Cauchy combination test)
GLMMMiRKAT(y = meta$y, X = cbind(meta$x1, meta$x2), id = meta$id, Ks = Ks, 
           model = "gaussian", method = "davies", omnibus = "perm")

# GLMM-MiRKAT with permutation test
GLMMMiRKAT(y = meta$y, X = cbind(meta$x1, meta$x2), id = meta$id, Ks = Ks, 
           model = "gaussian")

## -----------------------------------------------------------------------------
 # Import microbiome data with binomial traits
 data(bindata)

# Extract microbiome and meta information
# (Microbiome is the same as for Gaussian outcomes)
otu.tab <- bindata$bin.otu.tab
tree <- bindata$bin.tree
meta <- bindata$bin.meta
 
# Run GLMM-MiRKAT
GLMMMiRKAT(meta$y, X = cbind(meta$x1, meta$x2), id = meta$id, 
           Ks = Ks, model = "binomial")

## -----------------------------------------------------------------------------
# Import microbiome data with Poisson traits
data(poisdata)

# Extract microbiome and meta information
# (Microbiome is again the same)
otu.tab <- poisdata$pois.otu.tab
tree <- poisdata$pois.tree
meta <- poisdata$pois.meta

# Run GLMM-MiRKAT
GLMMMiRKAT(meta$y, X = cbind(meta$x1, meta$x2), id = meta$id, 
           Ks = Ks, model = "poisson")

