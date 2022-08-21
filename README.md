# Evaluation for Computational Methods

### Presented by Maria Camila Tavera Cifuentes

Welcome to the file that will guide you through this project! It was specially made for you and I hope you find it organized and useful.

First, here is presented the project structure.

## Project structure

    project/
    *    ├── .gitignore
    *    ├── data/
               └── raw
                    └── ClinicalTraits.csv
                    └── LiverFemale3600.csv
    *    ├── docs/
               └── project_report.html
               └── project_report.Rmd
               └── ref.bib
    *    ├── figs/
               └── modulesColorsDendogram.png
               └── phenoModRelation.png
               └── powerST.png
               └── sampleClusterDendogram.png
               └── traitHeatmap.png
    *    ├── Mouse_Liver_Gene_Exp.Rproj
    *    ├── output/
               └── femaleMouseTOM-block.1.RData
    *    ├── R/
             └── 01_read_data.R
             └── 02_packages_and_parameters.R
             └── 03_data_filtering.R
             └── 04_dendogram_heatmap.R
             └── 05_network_construction.R
             └── 06_pheno_mod_relation.R
             └── main.R

    *    ├── README.md

Now, let see what's inside the folders:

-   data/raw:

    Data comes from an study on a F2 mouse intercross, which objective was to find the relationship between gene expression patterns and related phenotype indicators in liver of 135 females.

    -   ClinicalTraits.csv: this is the phenotype data.

    -   LiverFemale3600.csv: this is the gene expression data.

-   figs:

    -   modulesColorsDendogram.png: this file is generated by the script *06_pheno_mod_relation.R.*

    -   phenoModRelation.png: this file is generated by the script *06_pheno_mod_relation.R.*

    -   powerST.png: this file is generated by the script *05_network_construction.R.*

    -   sampleClusterDendogram.png: this file is generated by the script *03_data_filtering.R.*

    -   traitHeatmap.png: this file is generated by the script *04_dendogram_heatmap.R.*

-   output:

    -   femaleMouseTOM-block.1.RData: this file is generated by the script *05_network_construction.R.*

-   R:

    -   01_read_data.R: this file reads the two .csv files located at *data/raw* folder.

    -   02_packages_and_parameters.R: this file loads all the libraries needed in order to run the codes. It is important to check this file before running any of the scripts because it has some instructions at the beginning.

    -   03_data_filtering.R: this file filters information for both gene expression and phenotype data sets.

    -   04_dendogram_heatmap.R: this file unifies the information found in the .csv files and presents a sample dendogram and trait heatmap.

    -   05_network_construction.R: this file executes an analysis of scale free topology for multiple soft thresholding powers, constructs a network and detects modules based on genes connectivity.

    -   06_pheno_mod_relation.R: this file creates a cluster dendogram and presents correlations between phenotypes and modules.

-   main.R: this file calls the other six scripts.

## Instructions

There are three steps to follow in order to run the codes and get the results that are expected.

1.  The first thing is to read the README. So, if you are here, you have almost completed the first step.
2.  The second step is (as it was said before) open the script *02_packages_and_parameters.R* and uncomment the three installation related lines, run these three and then, comment them back again.
3.  The third and last step is running the script *main.R*.

That's it! If you follow the instructions carefully, at the end of step 3 you should have all the figures and output presented on the project structure.

However, if something goes wrong, do not hesitate on contacting any person of the team.

Here is presented de contact information:

**Email:** mc.tavera\@uniandes.edu.co
