# Evaluation for Computational Methods ICTP-Serrapilheira 
# Script to create a coexpression network with WGCNA and get correlations with phenotypic features

#-----LIBRARIES INSTALLATION
#The installation lines should be uncomment:

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("DESeq2")
#BiocManager::install("WGCNA")
#BiocManager::install(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "ggplot2", "dplyr", "lattice"))

library(BiocManager)
library(impute)
library(GO.db)
library(preprocessCore)
library(org.Mm.eg.db)
library(DESeq2)
library(WGCNA)
library(ggplot2)
library(dplyr)
library(lattice)
library(matrixStats)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(survival)

#----PARAMETERS FOR PACKAGE WELL FUNCTIONING 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#----READING DATA

files_path <- list.files(path = "data/raw",
                         pattern = ".csv",
                         full.names = TRUE)

file_names <- gsub(".csv", "", basename(files_path), fixed = TRUE)
for (i in 1:length(files_path)) {
  data <- read.csv(files_path[[i]])
  assign(file_names[i], data)
}

#---EXPRESSION MATRIX HANDLING
#First, we create a subset data frame, taking into account expression data of
#each female mouse. Corresponding information goes from column 9 to the end.
#Then we transpose the information and the new columns names are the rows of
#substanceBXH original column.

LiverDF <- as.data.frame(t(LiverFemale3600[, -c(1:8)]))
names(LiverDF) <- LiverFemale3600$substanceBXH
rownames(LiverDF) <- names(LiverFemale3600)[-c(1:8)]

#---FILTERING
#WGCNA has a function  that checks for missing entries or low values, so the 
#next step is to use it and filter. It takes as it first parameter the expression matrix
#where samples (mice) are on the rows and on the columns are the genes. 
#The function checks if the samples and genes are considered good, and it return booleans.

gsg <- goodSamplesGenes(LiverDF, verbose = 3)

#Then we check if everything is good and in case it is not, we redefine the 
#dataframe so we just takes the good genes and samples

gsg$allOK

if (!gsg$allOK)
{
  LiverDF <- LiverDF[gsg$goodSamples, gsg$goodGenes]
}

#Now we check for outliers and we take them off using hierarchical clustering
#using "average" as the agglomeration method.




###MEJORAR?--------------------------TODO---------
sampleTree = fastcluster::hclust(dist(LiverDF), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red")


#After looking at the dendogram it is necessary to find the cut height, which
#is basically the outliers starting point

cutHeight <- 15 
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)

#cutreeStatic returns a vector with the cluster of each gene. In this vector it
#is seen that the only sample that is grouped in a different cluster than the 
#others F2_221, which was the outlier seen in the dendogram. in this order
#of ideas, we decide to stay with all samples but F2_221.

keepSamples <- (clust==1)
LiverDF_filtered <- LiverDF[keepSamples, ]
nGenes <- ncol(LiverDF_filtered)
nSamples <- nrow(LiverDF_filtered)

#NOW WE ARE GOING TO MANAGE PHENOTYPES' FILES 

phenotypes <- ClinicalTraits
names(phenotypes)

#On the dataset are columns associated to comments and notes, and they not 
#give relevant information, so they need to be deleted.
phenotypes_filtered <- phenotypes[, -c(31, 16)]

#Also we forget about other information associated columns, and we stay wiht
#phenotypes' associated columns and the column that has every mouse id.

phenotypes_filtered <- phenotypes_filtered[, c(2, 11:36)]
names(phenotypes_filtered)

#Now we organize data in a data frame analogous to the expression matrix
samplesExpMatrix <- rownames(LiverDF_filtered)


####--- TODO ------------------------- ACA SE DEBE PODER USAR ALGO DE LO DE TIDY

#We look which mice are in both data frames.
matchingSamples <- match(samplesExpMatrix, phenotypes_filtered$Mice)

#We create a new dataframe that has phenotype information just for the mice
#that are in both data frames.

matchingPhenotypeDF <- phenotypes_filtered[matchingSamples, -1]

#We change the rownames by mice IDs
rownames(matchingPhenotypeDF) <- phenotypes_filtered[matchingSamples, 1]

#We generate the dendogram again without outliers and showing phenotypes

sampleTree2 = hclust(dist(LiverDF_filtered), method = "average")

#We declare color codes for phenotypes: low value is white, high value is red 
#and nule is grey.


#---------------------TODO---mirar si se puede hacer más lindo

traitColors = numbers2colors(matchingPhenotypeDF, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(matchingPhenotypeDF), 
                    main = "Sample dendrogram and trait heatmap")

#Now we choose coexpression powers to calculate adjacency, as softsoft-tresholds,
#bases on network topology.

powers = c(c(1:10), seq(from = 12, to=24, by=2))
sft = pickSoftThreshold(LiverDF_filtered, dataIsExpr = TRUE, powerVector = powers,  corFnc = cor,networkType = "signed")
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") +
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  theme_classic()

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#PUNTO 2: Escoger un potenciador soft treshold para la produccion de la matriz de adyacencias y la red,
#basado en el umbral de R^2 establecido en la grafica
potenciador = 12 #cambiar de acuerdo a la grÃ¡fica

#Con este codigo se corre la red y la asignacion de modulos como clusters de genes 
#basados en su conectividad, todo se guarda como archivos de R aparte para manejar 
#la complejidad de espacio y tiempo de la computacion 
red = blockwiseModules(LiverDF_filtered, power = potenciador,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
#podemos observar cuantos modulos y cuantos genes por modulo quedaron en el proceso de clustering
table(red$colors)

#grafica que genera un dendrograma de los modulos con colores especificos para cada uno
sizeGrWindow(12, 9)
mergedColors = labels2colors(red$colors)
plotDendroAndColors(red$dendrograms[[1]], mergedColors[red$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


etiquetasModulos = red$colors
coloresModulos = labels2colors(red$colors)

#a data frame containing module eigengenes of the found modules (given by colors)

MEs = red$MEs
geneTree = red$dendrograms[[1]]

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(LiverDF_filtered, coloresModulos)$eigengenes

#los reorganiza para que los similares por correlacion esten uno al lado del otro

MEs = orderMEs(MEs0)


CorModuloFenotipo = cor(MEs, matchingPhenotypeDF, use = "p")
ValorPModuloFenotipo = corPvalueStudent(CorModuloFenotipo, nSamples)


#Esta grafica presenta las correlaciones y valores p entre phenotypes y clusters,
#lo cual permitira escoger los mas relevantes
textMatrix =  paste(signif(CorModuloFenotipo, 2), "\n(",
                    signif(ValorPModuloFenotipo, 1), ")", sep = "");
dim(textMatrix) = dim(CorModuloFenotipo)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = CorModuloFenotipo,
               xLabels = names(matchingPhenotypeDF),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Relaciones fenotipo-cluster"))

#PUNTO 3: Basado en la grÃ¡fica anterior escoger un modulo de interes diferente al default
#de acuerdo a su relacion con ciertos phenotypes y trabajar con este en adelante; ademas
#escoger un rasgo fenotipico con el cual se quiera trabajar, diferente a peso en este caso,
#cambiar el color para que quede ese modulo


modulo = "blue" #cambiar el color del modulo por el escogido
rasgo = as.data.frame(matchingPhenotypeDF$Glucose_Insulin) #en este caso se escoge el atributo weight_g, escoger otro
names(rasgo) = "rasgo" #no cambiar este String, el rasgo que escogio tendra este nombre en los atributos
nombresModulos = substring(names(MEs), 3)

pertenenciaGenModulo = as.data.frame(cor(LiverDF_filtered, MEs, use = "p"));
valorMMP = as.data.frame(corPvalueStudent(as.matrix(pertenenciaGenModulo), nSamples));

names(pertenenciaGenModulo) = paste("MM", nombresModulos, sep="");
names(valorMMP) = paste("p.MM", nombresModulos, sep="");

significanciaGenFenotipo = as.data.frame(cor(LiverDF_filtered, rasgo, use = "p"));
valorGSP = as.data.frame(corPvalueStudent(as.matrix(significanciaGenFenotipo), nSamples));

names(significanciaGenFenotipo) = paste("GS.", names(rasgo), sep="");
names(valorGSP) = paste("p.GS.", names(rasgo), sep="");


columna = match(modulo, nombresModulos);
genesModulo = coloresModulos==modulo;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(pertenenciaGenModulo[genesModulo, columna]),
                   abs(significanciaGenFenotipo[genesModulo, 1]),
                   xlab = paste("Module Membership in", modulo, "module"),
                   ylab = "Gene significance for glucose insulin",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulo)

names(LiverDF_filtered)[coloresModulos==modulo] #cambiar por el color de modulo escogido

#Se lee el archivo de anotaciones de genes
annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(LiverDF_filtered)
probes2annot = match(probes, annot$substanceBXH)
sum(is.na(probes2annot))
#esta suma debe dar 0 ya que son los genes sin anotacion

#Crear dataframe con informacion de anotacion y valores calculados
infoGenes0 = data.frame(substanceBXH = probes,
                        geneSymbol = annot$gene_symbol[probes2annot],
                        LocusLinkID = annot$LocusLinkID[probes2annot],
                        moduleColor = coloresModulos,
                        significanciaGenFenotipo,
                        valorGSP)
ordenModulos = order(-abs(cor(MEs, rasgo, use = "p")))
for (mod in 1:ncol(pertenenciaGenModulo))
{
  nombresAnteriores = names(infoGenes0)
  infoGenes0 = data.frame(infoGenes0, pertenenciaGenModulo[, ordenModulos[mod]], 
                          valorMMP[, ordenModulos[mod]])
  names(infoGenes0) = c(nombresAnteriores, paste("MM.", nombresModulos[ordenModulos[mod]], sep=""),
                        paste("p.MM.", nombresModulos[ordenModulos[mod]], sep=""))
}
ordenGenes = order(infoGenes0$moduleColor, -abs(infoGenes0$GS.rasgo)) #
infoGenes = infoGenes0[ordenGenes, ]

#se imprime un dataframe con la informacion de cada gen anotada, su pertenencia a cierto
#modulo y la informacion de correlacion gen fenotipo
write.csv(infoGenes, file = "infoGenes.csv")

#ahora se hara un enriquecimiento de terminos GO para encontrar las funciones de los genes
#del cluster y relacionarlas con los phenotypes relacionados

#Se realiza la anotacion GO para todos los genes
LLIDs = annot$LocusLinkID[probes2annot]
modLLIDs = LLIDs[genesModulo]
GOenr = GOenrichmentAnalysis(coloresModulos, LLIDs, organism = "mouse", nBestP = 10)

#This function is deprecated and will be removed in the near future. 
#We suggest using the replacement function enrichmentAnalysis 
#in R package anRichment, available from the following URL:
# https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/

tab = GOenr$bestPTerms[[4]]$enrichment

#se filtran los mejores resultados para el modulo escogido
tab_modulo = subset(tab, tab$module == modulo)
tab_modulo$termDefinition
#se imprime y se puede observar las funciones en la ultima columna
write.table(tab_modulo, file = "EnriquecimientoGOModulo.csv", sep = ",", quote = TRUE, row.names = FALSE)

#Exportar a Cytoscape para visualizar el modulo escogido
TOM = TOMsimilarityFromExpr(LiverDF_filtered, power = potenciador)
modulos = c(modulo)
probes = names(LiverDF_filtered)
enModulo = is.finite(match(coloresModulos, modulos))
modProbes = probes[enModulo]
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)]
modTOM = TOM[enModulo, enModulo]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-",
                                                paste(modulos, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-",
                                                paste(modulos, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = coloresModulos[enModulo])



