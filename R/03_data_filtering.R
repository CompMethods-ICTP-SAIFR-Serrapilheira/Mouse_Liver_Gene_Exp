# Evaluation for Computational Methods ICTP-Serrapilheira 

#This code filters information for both gene expression and phenotype datasets.

#------------------------EXPRESSION MATRIX (EM) HANDLING------------------------

# First, we create a subset data frame, taking into account expression data of
# each female mouse. Corresponding information goes from column 9 to the end.
# Then we transpose the information and the new columns names are the rows of
# substanceBXH original column.

LiverDF <- as.data.frame(t(LiverFemale3600[, -c(1:8)]))
names(LiverDF) <- LiverFemale3600$substanceBXH
rownames(LiverDF) <- names(LiverFemale3600)[-c(1:8)]

#---------------------------------EM FILTERING----------------------------------

# WGCNA has a function  that checks for missing entries or low values, so the
# next step is to use it and filter. It takes as it first parameter the expression
# matrix where samples (mice) are on the rows and on the columns are the genes. 
# The function checks if the samples and genes are considered good and it returns booleans.

gsg <- goodSamplesGenes(LiverDF, verbose = 3)

# Then we check if everything is good and in case it is not, we redefine the 
# dataframe so we just take the good genes and samples

gsg$allOK

if (!gsg$allOK){
  LiverDF <- LiverDF[gsg$goodSamples, gsg$goodGenes]
}

# Now we check for outliers and we take them off using hierarchical clustering
# using "average" as the agglomeration method.

sampleTree <- fastcluster::hclust(dist(LiverDF), method = "average")
png("figs/sampleClusterDendogram.png",width=1200, height=600)
par(cex = 0.6, mar = c(1,6,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", cex.lab = 1.2, 
     cex.axis = 1.2, cex.main = 1.5)
abline(h = 15, col = "red")
dev.off()

# After looking at the dendogram it is necessary to find the cut height, which
# is basically the outliers' starting point.

cutHeight <- 15 
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)

#cutreeStatic returns a vector with the cluster of each gene. In this vector it
#is seen that the only sample that is grouped in a different cluster (0) than the 
#others is F2_221, which was the outlier seen on the dendogram. In this order
#of ideas, we decided to stay with all samples but F2_221.

keepSamples <- (clust==1)
LiverDF_filtered <- LiverDF[keepSamples, ]
nGenes <- ncol(LiverDF_filtered)
nSamples <- nrow(LiverDF_filtered)

#---------------------------PHENOTYPES' FILE HANDLING---------------------------

phenotypes <- ClinicalTraits
names(phenotypes)

# On the dataset are columns not associated to phenotypic features (comments,
# notes, other mice info), so they are erased.

phenotypes_filtered <- phenotypes[, -c(31, 16, 1, 3:10)]
rownames(phenotypes_filtered) <- phenotypes_filtered[,1]
phenotypes_filtered <- phenotypes_filtered[,-1]