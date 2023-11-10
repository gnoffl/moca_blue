#rdf5_seqlet_pattern_v1.0.R
#DOCU:
#  This script automates the extraction and transformation of sequence patterns from an HDF5 file to a structured text format.
#Expression predictive motifs (EPMs) have distinct range from where they can be used effectively for gene expression levels.
#This scripts finds and and extracts the EPMs preferred sequence range, so that these can be used for optimized searches in genomes.
#It is designed for customizable naming conventions that can be used for EPM nomenclature and will be applied in downstream analyses of the moca_blue suite.
#The script generates a formatted output file in the specified directory.
#SETUP:
#  The script initializes the environment and sets various parameters:
#  NAME0: Name identifier for the output file, SPEC: Specifies the species (e.g., "Soly" for Solanum lycopersicum), MODEL: Specifies the model (e.g., "S0").
#FILE1: The name of the input HDF5 file.
#To execute and store the script, a "out" directory must be present in the parent directory for storing output files accessed by other scripts.
#FILEPATHs:
#  dirpath_in: The directory path where the input HDF5 file is located, dirpath_out: The directory path where the output sequences will be stored.
#LIBRARY DEPENDENCIES:
#  The script requires the rhdf5 and tidyr R packages. Uncomment the install.packages() line if the packages are not installed.
#PROCESSING STEPS:
#[1] - The script opens the input HDF5 file, reads metacluster information, and extracts sequence patterns for further processing.
#[2] - It iterates through metaclusters, retrieves sequence patterns, and organizes them into a structured format.
#[3] - The script formats and names the sequences based on a defined nomenclature and the number of patterns associated with each metacluster.
#[4] - Finally, the formatted sequences are written to a text file with a name constructed from specified parameters.

# 2023-11-10 Dr. Simon M. Zumkeller                                           # # #
                                                                    #   #    # # # #
########################## [SETUP] ############################## #  # # # # # #  # #
#install.packages("plyr")                                                  # # #
#library(plyr)
library(tidyr)
library(rhdf5)
###################################################################################
setwd("/home/ibg-4/Desktop/Rhome/moca_blue/mo_range")
###################################################################################
NAME0="rdf5_seqlet_pattern"
SPEC="Soly"
MODEL="S0"  # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
#######################################################
FILE1= "solanum_modisco.hdf5"
###################################################################################
dirpath_in = "../0MOTIFS/MODISCO_RAW/MODISCO_SSR"
dirpath_out = "./out"
# Define the pattern names to iterate over

###################################################
#motifs in metacluster 0 -automatize this step
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#X=20
#motifs in metacluster 1 -automatize this step
#Y=10
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    [1]
h5file <- H5Fopen(file.path(
  dirpath_in,
  FILE1), "H5F_ACC_RDONLY")
h5ls(h5file)
metacluster_group <- h5read(h5file, "metacluster_idx_to_submetacluster_results")
# loop through the metaclusters 0 and 1
for (i in c(0, 1)) {
  metacluster <- metacluster_group[[paste0("metacluster_", i)]]
}
#######################################################   [2]
# loop through the metaclusters 0 and 1
for (i in names(metacluster_group)) {
  metacluster <- metacluster_group[[i]]
  patterns = metacluster[['seqlets_to_patterns_result']][['patterns']]
}
# Define the pattern names to iterate over
########################################################
length(patterns[['all_pattern_names']])
X = length(metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]])-2
Y = length(metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]])-2
X1 = X+Y
Y1 = X*2+Y*2

##############################################             METACLUSTER_0    ###### [3]
# Initialize a list to store the results
ls_list <- list()
ls_list1 <- list()
seqlets_all_mc0 <- data.frame()
seqlets_all_mc1 <- data.frame()

for (i in 0:X) {
  pattern_name <- paste0("pattern_", i)
  seqletls <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  ls <- as.data.frame(seqletls)
  ls_list[[pattern_name]] <- ls
  seqlets_i <- as.data.frame(ls_list[[paste0("pattern_", i)]][["seqletls"]])
  colnames(seqlets_i) <- c("seqlets")
  seqlets_i$pattern <- paste0("pattern_", i)
  seqlets_all_mc0 <- rbind(seqlets_all_mc0, seqlets_i)
}
seqlets_all_mc0$metacluster <- c("metacluster_0")

##############################################        METACLUSTER_1  #############

for (i in 0:Y) {
  pattern_name <- paste0("pattern_", i)
  seqletls <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  ls <- as.data.frame(seqletls)
  ls_list[[pattern_name]] <- ls
  seqlets_i <- as.data.frame(ls_list[[paste0("pattern_", i)]][["seqletls"]])
  colnames(seqlets_i) <- c("seqlets")
  seqlets_i$pattern <- paste0("pattern_", i)
  seqlets_all_mc1 <- rbind(seqlets_all_mc1, seqlets_i)
}
seqlets_all_mc1$metacluster <- c("metacluster_1")

############################################################################ ###### 

seqlet_mc01 <- rbind.data.frame(seqlets_all_mc0, seqlets_all_mc1)

df <- seqlet_mc01 %>%
  mutate(example = NA, start = NA, end = NA, rc = NA) %>%
  separate(col = seqlets, into = c("example", "start", "end", "rc"), sep = "[,]")

df$example <- gsub("example:", "", df$example)
df$start <- gsub("start:", "", df$start)
df$end <- gsub("end:", "", df$end)
df$rc <- gsub("rc:", "", df$rc)

############################################################################ ###### [4] 
file_path_out <- file.path(dirpath_out, paste0(NAME0,SPEC,MODEL)) 

write.csv(df, file = paste0(file_path_out,".txt"), row.names = FALSE)
