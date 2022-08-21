# Evaluation for Computational Methods ICTP-Serrapilheira 

# This code creates a cluster dendogram and based on the module eigengenes and
# network colors is able to present correlations between phenotypes and modules.

# Now we graph a dendogram of the modules with specific colors for each of them.

modulesColors <- labels2colors(network$colors)
modulesLabels <- network$colors
png("figs/modulesColorsDendogram.png",width=1500, height=800)
plotDendroAndColors(network$dendrograms[[1]], modulesColors[network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Now we create a variable to save a data frame containing module eigengenes of
# the found modules (given by colors).

MEs <- network$MEs

# Now we recalculate module eigengenes with the assigned module colors.

MEs0 <- moduleEigengenes(LiverDF_filtered, modulesColors)$eigengenes

# Now reorganize so most similar ones (by correlation) stay close to each other.

MEs_Final <- orderMEs(MEs0)
corModPheno <- cor(MEs_Final, matchingPhenotypeDF, use = "p")
pValsCor <- corPvalueStudent(corModPheno, nSamples)

# Next graph presents correlations between phenotypes and modules
# which allows us to pick the most relevant ones.

textMatrix <- paste(signif(corModPheno, 2), "\n(",
                    signif(pValsCor, 1), ")", sep = "")

dim(textMatrix) <- dim(corModPheno)
png("figs/phenoModRelation.png",width=1200, height=500)
par(mar = c(10, 12, 3, 3))
labeledHeatmap(Matrix = corModPheno,
               xLabels = names(matchingPhenotypeDF),
               yLabels = names(MEs_Final),
               ySymbols = names(MEs_Final),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Phenotype-module relationships"))
dev.off()
