
### Workflow running the MTMM with the full sequencing data from Arabidopsis 1001 Genomes project 
### updated version 16.08.19
### This version uses ASREML-R versiion 4 

## first source the needed scripts and load libaries
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
#scriopts to source. All scripts can be found in the github folder 
source('scripts/emma.r')
source('scripts/mtmm_estimates_as4.r')
source('scripts/plots_gwas.r')
source('scripts/plot_mtmm.r')
source('scripts/mtmm_cluster.r')
source('scripts/mtmm_part2.r')
source('scripts/gwas.r')
## load your Phenotype  Y: a n by m matrix, where n=number of individuals and the first column contains the individual names (colname = 'ecotype_id') , the second the phenotypic values for the first trait and the third column the phenotypic values of the second trait
## load the Kinship matrix
## to test the script load the samle data 
#  load('/data/MTMM_SAMPLE_DATA.Rdata')

##run the mtmm_estimate script 
mtmm_estimates(Y,k=2,l=3,K,method='default',only.vca=FALSE) 

## this script will generate all the needed data for the GWAS using the estimates from the asreml call
## it will save all data as an .rda file named [Y1]_[Y2]_mtmm_estimates.rda
## you need to load this data before you can continue.
k=2
l=3
mydata<-paste(colnames(Y)[k],colnames(Y)[l],'mtmm_estimates.rda',sep='_')

load(mydata)
## create the name of your output file 
out.name=paste(colnames(Y)[k],'_',colnames(Y)[l],'_mtmm_results',sep='')

## mtmm_cluster is a wrapper script to run the MTMM with data from the Arabidopsis 1001 Genomes project in a sequential way
## you can get the sequencing data here https://go.uniwue.de/2029data

mtmm_cluster(out.name,incl.singleGWAS=T)

## this script will generate Rdata with the pvalues of the full test, the trait specific test and the trait common test, additionally the single GWAS analyses can be included
 
## for other data then the ARabidopsis full-sequence data just run the mtmm_part2 script with your own genotype data (X)
# results<-mtmm_part2(X)
