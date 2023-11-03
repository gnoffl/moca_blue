# rdf5_get_cwms_per_pattern_v1.0.R
# DOCU:
# This script is designed to automate the extraction and transformation of motifs from HDF5 files to JASPAR format. 
# It allows for customization of naming conventions and can be used in various bioinformatics and motif analysis tasks.
# The script will generate JASPAR-formatted motif files in the specified output directory.
# Below is the more detailed documentation for scripts modules:
#
# SETUP:
# The script begins by setting up the environment and specifying various parameters:
# NAME0: Name identifier for source file, SPEC: Specifies the species (e.g., "Zema" for Zea mays), MODEL: Specifies the model (e.g., "S0", "S1" or "S0+M0")
# FILE1: The name of the input HDF5 file.
# To run and store the script the parent directory requires a daughter directory named "out" where the output is stored and other scripts of the moca_blue suite will access it.
# FILEPATHs: 
# dirpath_in: The directory path where the input HDF5 file is located, dirpath_out: The directory path where the output motifs will be stored.
# LIBARY DEPENDENCIES: 
# The script loads several R packages, including rhdf5, and tidyr. 
# 
# PROCESSING STEPS:
# [1] - The script reads the input HDF5 file and extracts metacluster information and patterns for further processing.
# [2] - The patterns extracted from the HDF5 file are transformed into CWM (Contribution Weight Matrices) format. 
#       This includes extracting forward and reverse matrices for each pattern. The contribution scores for each pattern are converted into weight matrices (CWM) 
#       by scaling them based on the number of seqlets associated with the pattern.
# [3] - The script assigns names to the motifs based on the proposed nomenclature and the number of seqlets associated with each pattern. 
#       It also ensures that the names are unique. Thourghout the whole DeepCRE analyses the given nomenclature is beeing used. For a detailed description please refer
#       to our corresponding manuscript: https://doi.org/10.21203/rs.3.rs-2873437/v1
# [4] - Finally, the transformed motifs in JASPAR format are written to an output file with a name constructed from the specified parameters.

# 2023-10-05 by Dr. Simon M. Zumkeller

########################## [SETUP] ##############################
setwd("/home/ibg-4/Desktop/Rhome/moca_blue/mo_nom")
###################### 
NAME0="rdf5_1"
SPEC="Zema"
MODEL="S0"
#######################################################
FILE1= "zea_modisco.hdf5"
#######################################################
dirpath_in = "../0MOTIFS/MODISCO_RAW/MODISCO_SSR"
dirpath_out = "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(rhdf5)
library(tidyr)
########################## [PROCESSING] ##############################
#                                                                           [1]
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
h5file <- H5Fopen(file.path(
  dirpath_in,
  FILE1), "H5F_ACC_RDONLY")
metacluster_group <- h5read(h5file,
                            "metacluster_idx_to_submetacluster_results")
#######################################################
for (i in names(metacluster_group)) {
  metacluster <- metacluster_group[[i]]
  patterns = metacluster[['seqlets_to_patterns_result']][['patterns']]
}
########################################################
length(patterns[['all_pattern_names']])
x0 = length(metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]])-1
x1 = length(metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]])-1
y01 = x0+x1
y02 = x0*2+x1*2
###################################################
pattern_names <- paste0("pattern_", 0:y01) ############################################## !!! MANUAL ADJ REQUIRED
matricesF0 <- list()
matricesF1 <- list()
matricesR0 <- list()
matricesR1 <- list()
for (pattern_name in pattern_names) {
  matrixF0 <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["fwd"]]
  matrixF1 <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["fwd"]]
  matrixR0 <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["rev"]]
  matrixR1 <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["task0_contrib_scores"]][["rev"]]
  matricesF0[[pattern_name]] <- matrixF0
  matricesF1[[pattern_name]] <- matrixF1
  matricesR0[[pattern_name]] <- matrixR0
  matricesR1[[pattern_name]] <- matrixR1
}
##################################################################################### PFM to PWM to CWM
#                                                                           [2]
seqletls_lengths_p1 <- list()
for (pattern_name in pattern_names) {
  seqletls <- metacluster_group[["metacluster_1"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  seqletls_lengths_p1[[pattern_name]] <- length(seqletls)
}
seqletls_lengths_p0 <- list()
for (pattern_name in pattern_names) {
  seqletls <- metacluster_group[["metacluster_0"]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
  seqletls_lengths_p0[[pattern_name]] <- length(seqletls)
}
  ###############                ASSIGN NOMENCLATURE            ###################
 ###############                ASSIGN NOMENCLATURE            ###################
###############                ASSIGN NOMENCLATURE            ###################
#                                                                           [3]
names(matricesF0) <- paste0(names(matricesF0),
                            "_p0m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesF0), 9))))
names(matricesF0) <- substr(names(matricesF0), nchar(names(matricesF0)) - 4, nchar(names(matricesF0)))
names(matricesF0) <- paste0(names(matricesF0),
                            "F")
names(matricesF0) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesF0))
###############                  ####################            ###################
names(matricesF1) <- paste0(names(matricesF1),
                            "_p1m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesF1), 9))))
names(matricesF1) <- substr(names(matricesF1), nchar(names(matricesF1)) - 4, nchar(names(matricesF1)))
names(matricesF1) <- paste0(names(matricesF1),
                            "F")
names(matricesF1) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesF1))
###############                  ####################            ###################
names(matricesR1) <- paste0(names(matricesR1),
                            "_p1m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesR1), 9))))
names(matricesR1) <- substr(names(matricesR1), nchar(names(matricesR1)) - 4, nchar(names(matricesR1)))
names(matricesR1) <- paste0(names(matricesR1),
                            "R")
names(matricesR1) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesR1))
###############                  ####################            ###################
names(matricesR0) <- paste0(names(matricesR0),
                            "_p0m",
                            sprintf("%02d",
                                    as.numeric(substring(names(matricesR0), 9))))
names(matricesR0) <- substr(names(matricesR0), nchar(names(matricesR0)) - 4, nchar(names(matricesR0)))
names(matricesR0) <- paste0(names(matricesR0),
                            "R")
names(matricesR0) <- paste0("epm_", SPEC, "_", MODEL, "_", names(matricesR0))
###############                  ####################            ###################   # ADD NUMBER OF SEQLETS TO NAMES !!!!
seqlets_count_p0 <- head(seqletls_lengths_p0, x0)

for (i in seq_along(seqlets_count_p0)) {
  name <- paste0(names(matricesF0)[i], "_", seqlets_count_p0[[i]])
  names(matricesF0)[i] <- name
}

for (i in seq_along(seqlets_count_p0)) {
  name <- paste0(names(matricesR0)[i], "_", seqlets_count_p0[[i]])
  names(matricesR0)[i] <- name
}

seqlets_count_p1 <- head(seqletls_lengths_p1, x1)

for (i in seq_along(seqlets_count_p1)) {
  name <- paste0(names(matricesF1)[i], "_", seqlets_count_p1[[i]])
  names(matricesF1)[i] <- name
}

for (i in seq_along(seqlets_count_p1)) {
  name <- paste0(names(matricesR1)[i], "_", seqlets_count_p1[[i]])
  names(matricesR1)[i] <- name
}
#################################################################################### contribution scores, to weitgh matrix
motifs <- c(matricesF0,matricesF1,matricesR0,matricesR1)
####################################################################################
for (i in seq_along(motifs)) {
  m0 <- motifs[i]
  m1 <- motifs[[i]]
  name <- names(m0)[1]
  seq_count <- sub(".*_([0-9]+)$", "\\1", name)
  nfcwm <- abs(m1)
  nfcwm <- round(as.numeric(seq_count)*(nfcwm/max(nfcwm)))
  motifs[[i]] <- abs(nfcwm)
}

####################################################################################
pfms<- array(unlist(motifs),dim = c(4, 14, y02))    # make correction here - not pfm but cwm!!!!! extracting from sequence gives the PPM ! position probability matrix!!!!
ls_pfms<- list()
for (idx in seq(1:y02)){
  ls_pfms[[idx]] <-pfms[, , idx]
}
ls_pfms_str <- lapply(ls_pfms, function(x) {
  apply(x, c(1, 2), as.character)
})
#####################################################################################
create_text <- function(m){
  res <- paste0(">motif", "\n")
  rows <- c('A', 'C', 'G', 'T')
  for(i in 1:nrow(m)){
    res <- paste0(res, paste0(rows[i],' ', paste0('[', paste(as.character(m[i, 1:ncol(m)]), collapse = "\t"), ']', "\n")))
  }
  return(res)
}
##########################                            ###############################
                          ############################
text <- lapply(ls_pfms, create_text)
for (idx in seq(1:y02)){
  text[[idx]] <- gsub("motif", paste0(names(motifs)[idx]), text[[idx]])
}
#############
#                                                                           [4]
file_path_out <- file.path(dirpath_out, paste0(NAME0,SPEC,MODEL))

writeLines(unlist(text), paste0(file_path_out,"_cwm-motifs.jaspar"))
#####################################################################################
