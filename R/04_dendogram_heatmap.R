# Evaluation for Computational Methods ICTP-Serrapilheira 

# This code unifies the information found in the .csv files and presents
# a sample dendogram and trait heatmap.

#-------------------------------MATCHING MATRICES-------------------------------

# Now we organize data in a data frame analogous to the expression matrix.

samplesExpMatrix <- rownames(LiverDF_filtered)
samplesPhenotype <- rownames(phenotypes_filtered)

# We look for the mice of the expression matrix dataset that are also on the 
# phenotype dataset.

matchingSamples <- match(samplesExpMatrix, samplesPhenotype)

# We create a new dataframe that has phenotype information just for the mice
# that are in both data frames.

matchingPhenotypeDF <- phenotypes_filtered[matchingSamples,]

#-----------------------SAMPLE DENDOGRAM AND TRAIT HEATMAP----------------------

# We generate the dendogram again without outliers and showing phenotypes

sampleTree2 <- hclust(dist(LiverDF_filtered), method = "average")

# We declare color codes for phenotypes: low value is white, high value is red 
# and nule is grey.

traitColors <- numbers2colors(matchingPhenotypeDF, signed = FALSE)
png("figs/traitHeatmap.png",width=1300, height=600)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(matchingPhenotypeDF), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
