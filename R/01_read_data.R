# Evaluation for Computational Methods ICTP-Serrapilheira 


# This code reads the .csv files that are in data folder. 
# The two files that will be read are LiverFemale3600.csv and ClinicalTraits.csv
# The first one presents expression data for female mice and the second one 
# presents phenotype related data.

files_path <- list.files(path = "data/raw",
                         pattern = ".csv",
                         full.names = TRUE)

file_names <- gsub(".csv", "", basename(files_path), fixed = TRUE)
for (i in 1:length(files_path)) {
  data <- read.csv(files_path[[i]])
  assign(file_names[i], data)
}

