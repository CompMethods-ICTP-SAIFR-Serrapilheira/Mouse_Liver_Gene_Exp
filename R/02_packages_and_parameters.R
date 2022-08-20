# Evaluation for Computational Methods ICTP-Serrapilheira 


# This code loads the packages that will be needed and establishes some 
# parameters for the future functions that will be used.

# The following  3 commented lines should be uncommented the first time the script is
# used. After that, they should be commented again.

# BiocManager::install("DESeq2")
# BiocManager::install("WGCNA")
# install.packages(c("foreach", "doParallel", "fastcluster", "dynamicTreeCut", "ggplot2", "dplyr", "lattice"))
library(BiocManager)
library(lattice)
library(foreach)
library(doParallel)
library(fastcluster)
library(WGCNA)
library(dynamicTreeCut)
library(dplyr)
library(ggplot2)

#----PARAMETERS FOR PACKAGE WELL FUNCTIONING 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()



