######################
library(rhdf5)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
######################################################
setwd("~/Desktop/Rhome/moca_blue/0MOTIFS")
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_"
SPEC="Arth"
MODEL="Dx" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
DATE= "20240205"
#######################################################
FILEs=c("Atha-0h_scores.h5",
        "Atha-1-5h_scores.h5",
        "Atha-3h_scores.h5",
        "Atha-6h_scores.h5",
        "Atha-12h_scores.h5",
        "Atha-24h_scores.h5"
        )
#######################################################
dirpath_in1 = "./"
dirpath_out = "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path1 = file.path(dirpath_in1,FILEs)
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################
# Define a function to process a file and store results
process_file <- function(file_name) {
  file_path <- file.path(dirpath_in1, file_name)
  
  h5file <- H5Fopen(file_path, "H5F_ACC_RDONLY")
  saliency_scores <- h5read(h5file, "contrib_scores")
  
  num_arrays <- dim(saliency_scores)[3]
  num_columns <- dim(saliency_scores)[2]
  result_matrix <- matrix(0, nrow = num_arrays, ncol = num_columns)
  
  for (i in 1:num_arrays) {
    current_array <- saliency_scores[,,i]
    column_sums <- colSums(current_array)
    result_matrix[i,] <- column_sums
  }
  
  # Store the result matrix in an object
  assign(gsub(".h5", "_result", file_name), result_matrix, envir = .GlobalEnv)
}

# List of files to process
files_to_process <- FILEs

# Process each file and store results
for (file_name in files_to_process) {
  process_file(file_name)
}

colav_0h <- colMeans(`Atha-0h_scores_result`)
colav_1_5h <- colMeans(`Atha-1-5h_scores_result`)
colav_3h <- colMeans(`Atha-3h_scores_result`)
colav_6h <- colMeans(`Atha-6h_scores_result`)
colav_12h <- colMeans(`Atha-12h_scores_result`)
colav_24h <- colMeans(`Atha-24h_scores_result`)
#######################################################
plot(colav_0h, type="l",
     col=rgb(1,0,0, alpha = 0.5),
     xlab="gene flanking regions",
     ylab="importance")

#rgb(239, 236, 236)
# Adding other lines to the same plot
lines(colav_1_5h, col=rgb(0.218,0.165,0.7, alpha = 0.3))
lines(colav_3h, col=rgb(0.142,0.199,0.210, alpha = 0.3))
lines(colav_6h, col=rgb(0.218,0.7,0.59, alpha = 0.5))
lines(colav_12h, col=rgb(0.7,0.71,0.90, alpha = 0.3))
lines(colav_24h, col=rgb(0,0.25,0.40, alpha = 0.3))
abline(h=0, col="black", lty=3)

legend("topright", inset = 0, xpd = TRUE, bty = "n",
       legend=c("0h","1.5h", "3h", "6h", "12h", "24h"), 
       col=c(rgb(1, 0, 0, alpha = 0.5),
             rgb(0.218, 0.165, 0.7, alpha = 0.7),
             rgb(0.142, 0.199, 0.210, alpha = 0.7),
             rgb(0.218,0.7,0.59, alpha = 0.7),
             rgb(0.7, 0.71, 0.90, alpha = 0.7),
             rgb(0, 0.25, 0.40, alpha = 0.7)), lty=1)