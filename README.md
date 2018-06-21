# MTMM - A mixed-model approach for genome-wide association studies of correlated traits in structured populations

##  Introduction

The MTMM function as published in [Nature Genetics](http://www.nature.com/ng/journal/v44/n9/full/ng.2376.html) currently don't support estimates on missing data and replicates. 
This is work in progress and will be accordingly updated here.

For questions and comments feel free to contact me: arthur.korte@gmail.com
 

## How to use

# To only perform a Variance Coponent Analysis use the mtmm_estimate.r script with the flag only.vca=T set
 
# Load libraries and source needed functions
# The AsREML package needs a valid license that can be obtained at  http://www.vsni.co.uk/software/asreml

library(lattice)
library(asreml)

# msm and nadiv librarys are used to estimate SE of the correlation estimates, only used if run=FALSE 
#library(msm)
#library(nadiv)

source('mtmm_function.r')
source('emma.r')

# load your data (Phenotype(Y),Genotype(X) and Kinship(K))
# note you can calculate K using the emma package K<-emma.kinship(t(X)), make sure to set colnames(K)=rownames(K)=rownames(X)

# alternativley load the sample data
load('data/MTMM_SAMPLE_DATA.Rdata')

# different options include method(default or errorcorrelation, include.single.analysis, calculate.effect.size (if TRUE, ###analysis is more time consuming) default for X is binary coding of 0 and 1, if your data are code  0,1 and 2 use ###gen.data='heterozygot',  run=FALSE will not perform the GWAS, but only output the correlation estimates (fast)
mtmm(Y,X,K,method='default',include.single.analysis=T,calculate.effect.size=T,gen.data='binary',exclude=T,run=T)

# the function outputs a list called results  ($phenotype ,$pvals, $statistics, $kinship)
output<-results$pvals

# manhattan plots
# default plots for include.single.analysis=T
par(mfrow=c(5,1),mar=c(3, 4, 1, 4))
plot_gwas(output,h=8)
plot_gwas(output,h=9)
plot_gwas(output,h=10)
plot_gwas(output,h=11)
plot_gwas(output,h=12)

#qq plots
par(mfrow=c(1,1),mar=c(3, 4, 1, 4))
qq_plot_all(output)
```
## Poster
* [Complextraits 2012](posters/poster_complextraits_2012_AK.pdf)
* [ICAR 2012](posters/poster_ICAR_2012_AK.pptx)

