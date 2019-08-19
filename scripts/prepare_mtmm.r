## dependencies of the MTMM script 
## ASREML library needs a valid license 
library(lattice)
library(asreml)
library(msm)
library(nadiv)
## libraries for single GWAS
library(foreach)
library(iterators)
library(parallel)
# libraries for plotting
library(ggplot2)
library(dplyr)  
#scriopts to source. All scripts can be found in the github folder. Be sure to set the right working directory
source('scripts/emma.r')
source('scripts/mtmm_estimates_as4.r')
source('scripts/plots_gwas.r')
source('scripts/plot_mtmm.r')
source('scripts/mtmm_cluster.r')
source('scripts/mtmm_part2.r')
source('scripts/gwas.r')