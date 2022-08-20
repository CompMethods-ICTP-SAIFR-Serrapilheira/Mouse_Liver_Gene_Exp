# Evaluation for Computational Methods ICTP-Serrapilheira 

# This code executes an analysis of scale free topology for multiple soft 
# thresholding powers in order to find the right power to create the network
# and also detect modules based on genes connectivity.

#------------------------------SOFT THRESHOLDING--------------------------------

# Now we need to pick an appropriate soft-thresholding power for network construction.

powers <- c(c(1:10), seq(from = 12, to=24, by=2))
sft <- pickSoftThreshold(LiverDF_filtered, dataIsExpr = TRUE, 
                         powerVector = powers,  corFnc = cor,
                         networkType = "signed")

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
ggsave("figs/powerST.png") 

# Based on the graph we decided to choose the soft-thresholding power = 12.

STpower <- 12

#------------------------------NETWORK CONSTRUCTION-----------------------------

# The next step is to construct the network and detect modules based on 
# genes connectivity (clusters). It is saved as an .RData file in output folder.

network <- blockwiseModules(LiverDF_filtered, power = STpower,
                            TOMType = "unsigned", minModuleSize = 30,
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "output/femaleMouseTOM",
                            verbose = 3)

# We want to see how many modules and genes per module remained after clustering process.
table(network$colors)
