#mo_cluster
# DOCU:
# This script is used to cluster motifs stored a weight matrices based on the Smith Waterman, Sandelin Wassermann (2004) or euclidian distance methods 
# and results are obtained as Newick tree (NWK tree) representation of the clustered motifs. 
# Additionally, it identifies sets of motifs with high similarity and separates them into two categories: 
# Motifs with highly similar counterparts and motifs without highly similar counterparts. 
# Below is the more detailed documentation for scripts modules:
#
# SETUP:
# The script starts by setting up some variables and libraries required for motif clustering.
# NAME0: Name identifier for source file, SPEC: Specifies the species (e.g., "Zema"), MODEL: Specifies the model (e.g., "S0", "S1" or "S0+M0"), TYPE: Specifies the type of motif file (e.g., "_cwm-motifs.jaspar")
# To run and store the script the parent directory requires a daugther directory named "out" where the output is stored and other scripts of the moca_blue suite will access it.
# FILEPATHs: 
# jaspar_file: Combines the setup variables to create the input motif file path, dirpath_in: Specifies the input directory path, dirpath_out: Specifies the output directory path.
# LIBARY DEPENDENCIES: 
# The script loads several R libraries including grid, TFBSTools, motifStack, universalmotif, ape, and ggtree.
# 
# PROCESSING STEPS:
# [1] - Motifs are read from the JASPAR motif file specified in jaspar_file using the read_jaspar function. The motifs are stored in the cwm1 list.
# [2] - Motif names are modified to include additional information. The script extracts the number of sites (nsites or consensus) from the motif name and appends it to the motif name.
# [3] - The motifs are converted into different formats, such as position weight matrices (PWM) and position count matrices (PCM) for later analyses depending on the clustering algorithm.
# [4] - The script summarizes the motifs in PCM format and writes the summary to a CSV file.
# [5] - Motifs are clustered using a selected method using the universalmotif library. Motifs can be compared using e.g. Sandelin Wassermann algorithm ("SW") or others. (See documentation https://bioconductor.org/packages/release/bioc/html/universalmotif.html)
# [6] - Highly similar motifs are identified by comparing the clustering results. Motifs that have highly similar counterparts are separated from those that don't.
# [6.1] The script also performs clustering using the Smith-Waterman method and saves the resulting NWK tree.
# [7] - From the clustering results, Highly similar (using the upper 5 percent quantile) motifs are identified . Motifs that have highly similar counterparts are separated from those that don't.
#       The comparison matrix, used here, is converted into a data frame. This data frame contains the motif similarity scores, with motifs as rows and columns. 
#       The quantile function is used to calculate the 5th percentile of the similarity scores in the data frame. 
#       This value represents a lower threshold for considering motifs as highly similar.
# [8] - Various outputs are written to files, including the NWK tree of clustered motifs, lists of motifs with and without highly similar counterparts, and a matrix of motif comparisons.

# 2023-10-05 by Dr. Simon M. Zumkeller
######################################################
################## [SETUP] ###########################
######################################################
#jaspar_file = paste0(NAME0,SPEC,MODEL,TYPE)
jaspar_file = "rdf5_ArthD0D6_cwm-motifs.jaspar"
#######################################################
string_to_remove1 <- "Arth_D0" #Define what EPMs should be compared by PCC similarity from here
string_to_remove2 <- "Arth_D6" # to here...
#######################################################
#######################################################
dirpath_out = "./out/"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## [COMMAND LINE ARGS] ##############################
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 3) {
  stop("Cannot provide more than 3 arguments! Usage: <script> [<working_directory>] [<jaspar_file>] [<dirpath_out>]", call.=FALSE)
}
if (length(args)>=3) {
  dirpath_out = args[3]
}
if (length(args)>=2) {
  jaspar_file = args[2]
}
if (length(args)>=1) {
  working_directory = args[1]
}
# print all arguments
cat("working_directory:", working_directory, "\n")
cat("jaspar_file:", jaspar_file, "\n")
cat("dirpath_out:", dirpath_out, "\n")


setwd(working_directory)

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to check and install packages
install_if_missing <- function(packages) {
  # Get list of installed packages
  installed_packages <- rownames(installed.packages())
  
  # Check which packages are missing
  missing_packages <- packages[!packages %in% installed_packages]
  
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    
    # Separate CRAN and Bioconductor packages
    bioc_packages <- c("TFBSTools", "motifStack", "universalmotif", "ggtree")
    cran_packages <- missing_packages[!missing_packages %in% bioc_packages]
    bioc_missing <- missing_packages[missing_packages %in% bioc_packages]
    
    # Install CRAN packages
    if (length(cran_packages) > 0) {
      install.packages(cran_packages, repos = "https://cran.r-project.org")
    }
    
    # Install Bioconductor packages
    if (length(bioc_missing) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.r-project.org")
      }
      BiocManager::install(bioc_missing)
    }
  }
}

# List of required packages
required_packages <- c("grid", "TFBSTools", "motifStack", "universalmotif", 
                      "ape", "ggtree", "ggplot2", "ggseqlogo", "Cairo")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.r-project.org")
      }
BiocManager::install(version = "3.21", ask=FALSE) # or omit 'version' for latest
# Install missing packages
install_if_missing(required_packages)

# Load libraries
library(grid)
library(TFBSTools)
library(motifStack)
library(universalmotif)
library(ape)
library(ggtree)
library(ggplot2)
library(ggseqlogo)
library(Cairo)
# ##################################################
cwm1 <- read_jaspar(jaspar_file)
output_path_base = file.path(dirpath_out, paste0(basename(jaspar_file)))
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  # Extract the motif name and the number after the last "_" underscore
  motif_name <- attr(cwm1[[i]], "name")
  nsites <- as.numeric(sub(".+_(\\d+)$", "\\1", motif_name))
  # Assign the nsites value to the "nsites" field of the motif object
  cwm1[[i]]["nsites"] <- nsites
}
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  # Extract the Total IC and the Consensus values from the motif object
  total_ic <- attr(cwm1[[i]], "icscore")
  total_ic_rounded <- round(total_ic, 1)
  consensus <- attr(cwm1[[i]], "consensus")
  # Combine the Total IC and the Consensus values separated by "_" to the motif name
  motif_name <- attr(cwm1[[i]], "name")
  new_motif_name <- paste0(motif_name, "_", total_ic_rounded, "_", consensus)
  # Assign the new motif name to the motif object
  attr(cwm1[[i]], "name") <- new_motif_name
}
##################################################
pwm_uni0<-convert_motifs(
  cwm1, class = "TFBSTools-PWMatrix")
pcm<-convert_motifs(
  cwm1, class = "motifStack-pcm")
##################################################
#ggseqlogo(pcm[[1]]@mat)
#pcm[[1]]@name
# Your existing code to generate the sequence logo
#seq_logo <- ggseqlogo(pcm[[1]]@mat)
#seq_logo <- seq_logo + theme_minimal() + theme(
#  plot.background = element_rect(fill = "white"),
#  panel.background = element_rect(fill = "white"),
#  legend.background = element_rect(fill = "white")
# Save the plot with custom width and height
#ggsave("your_plot.png", plot = seq_logo, width = 4831, height = 592, units = "px")
pcm_length <- length(pcm)  # Getting the length of the pcm object
for (i in 1:pcm_length) {
  seq_logo <- ggseqlogo(pcm[[i]]@mat)
  seq_logo <- seq_logo + theme_minimal() + theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white")
  )
  # Save the plot with custom width and height
  file_name <- paste0(pcm[[i]]@name, ".png")
  file_path <- file.path(dirpath_out, file_name)
  ggsave(file_path, plot = seq_logo, width = 4831, height = 592, units = "px")
}

#ignoring the following lines, as i just want to visualize the motifs

# ##################################################
# sum<-as.data.frame(summarise_motifs(pcm))
# write.csv(sum, file = paste0(output_path_base,"summary.txt"))
# ##################################################
# ##########################################################
# c<-compare_motifs(cwm1, method = "SW")
# c0<-as.data.frame(c)
# c1 <- as.matrix(c0)
# lower_percentile <- quantile(c1, probs = 0.05)
# #lower_percentile <- max(c1)*0.1
# upper_percentile <- quantile(c1, probs = 0.95)
# #upper_percentile <- max(c1)*0.9
# ###
# cX <- c0[!grepl(string_to_remove1, rownames(c0))]
# cX1 <- as.data.frame(t(cX))
# cX2 <- cX1[!grepl(string_to_remove2, rownames(c0))]
# ##############         ####################      ##################
# cY <- c0[!grepl(string_to_remove2, rownames(c0))]
# cY1 <- as.data.frame(t(cY))
# cY2 <- cY1[!grepl(string_to_remove1, rownames(c0))]
# ##############         ####################      ##################
# cXX2 <- cX2-upper_percentile
# cXX3 <- as.data.frame(ifelse(cXX2 < 0, 0, 1))
# cXX3$sum <- rowSums(cXX3)

# cYY2 <- cY2-upper_percentile
# cYY3 <- as.data.frame(ifelse(cYY2 < 0, 0, 1))
# cYY3$sum <- rowSums(cYY3)
# epms_without_highly_similar_counterpartsX <- rownames(cXX3[cXX3$sum < 1, ])
# epms_with_highly_similar_counterpartsX <- rownames(cXX3[cXX3$sum >= 1, ])
# epms_without_highly_similar_counterpartsY <- rownames(cYY3[cYY3$sum < 1, ])
# epms_with_highly_similar_counterpartsY <- rownames(cYY3[cYY3$sum >= 1, ])
# epms_without_highly_similar_counterparts <- c(epms_without_highly_similar_counterpartsX,
#                                                   epms_without_highly_similar_counterpartsY)
# epms_with_highly_similar_counterparts <- c(epms_with_highly_similar_counterpartsX,
#                                               epms_with_highly_similar_counterpartsY)
# #########################################################################
# write.table(epms_without_highly_similar_counterparts, file = paste0(output_path_base, "_epms_without_highly_similar_counterparts-SW.csv"), sep = "\t", col.names = NA, quote = FALSE)
# write.table(epms_with_highly_similar_counterparts, file = paste0(output_path_base, "_epms_with_highly_similar_counterparts-SW.csv"), sep = "\t", col.names = NA, quote = FALSE)
# write.table(c0, file = paste0(output_path_base, "_matrix-SW.csv"), sep = "\t", col.names = NA, quote = FALSE)
# #assign scores from data.frame to branches for selection
# comp_1 <- 1-c
# comp_1 <- as.dist(comp_1)
# #labels <- attr(comp_1, "Labels")
# comp_2 <- ape::as.phylo(hclust(comp_1))
# comp_2[["edge.length"]] <- comp_2[["edge.length"]]+1
# comp_2[["edge.length"]]
# # Create a rooted phylo object
# phylo_tree <- as.phylo(comp_2)
# # Save the tree with positive edge lengths
# write.tree(comp_2, file = paste0(output_path_base, "-Sandelin-Wassermann.nwk"))
# #################################################
# c_pcm<-clusterMotifs(pcm, method = "Smith-Waterman") ### !!! TIME TO GET A COFFEEE !!! ###
# hc<- c_pcm
# motifs<-pcm[hc$order]
# ##################################################
# write.tree(as.phylo(c_pcm), file = paste0(output_path_base, "-Smith-Waterman.nwk"))

# #par(cex = 0.2)##############################################################
# #pdf("output3.pdf", width = 30, height = 40)
# #plotMotifLogoStackWithTree(pcm, c_pcm,
# #                           treewidth = 1/8,
# #                           trueDist = TRUE)
# #dev.off()
# # Get the number of elements in pcm
# num_elements <- length(pcm)
# # Initialize an empty list to store filtered elements
# pcm_filtered <- vector("list", length = num_elements)
# # Loop through each element in pcm and filter
# for (i in 1:num_elements) {
#   if (!grepl("R_", pcm[[i]]@name)) {
#     pcm_filtered[[i]] <- pcm[[i]]
#   }
# }
# # Remove NULL elements (objects where @name contains "R_")
# pcm_filtered <- pcm_filtered[!sapply(pcm_filtered, is.null)]
# browseMotifs(
#   pcm_filtered,
#   layout = "cluster",
#   nodeRadius = 5,
#   baseHeight = 30,
#   baseWidth = 12,
#   xaxis = TRUE,
#   width = 100,
#   height = 1700
# )
# ##############################################################

# #pcm2 <- pcm
# #pcm2<- pcm2[!sapply(pcm2, function(x) grepl("R_", x$name))]
# #sapply(pcm2, class)
# #for (i in seq_along(pcm2)) {
# #  pcm2[[i]]$name <- substr(pcm2[[i]]$name, 1, 23)
# #}
# #str(pcm2)
# #pcm_hcNW <- clusterMotifs(pcm2, method = "Needleman-Wunsch")
# #phylo_tree <- ade4::hclust2phylog(pcm_hcNW)
# #phylo_tree <- as.phylo(pcm_hcNW)
# #plot(as.dendrogram((phylo_tree)))
# #plotMotifStackWithPhylog(phylo_tree, pcm2)

# #                               f.phylog = 0.3,
#                              #  cleaves = 1,
#  #                               cnodes = 0,
# #                               labels.leaves = names(phylo_tree$leaves),
# #                               clabel.leaves = 1,
# #                               labels.nodes = names(phylo_tree$nodes),
# #                               clabel.nodes = 0,
# #                               font = "Helvetica-Bold",
# #                             #  ic.scale = TRUE
# #)

# #plot(plot)
# #plotMotifLogoStackWithTree(pcm2,
# #                         pcm_hcNW,
# #                         treewidth= 1/4,
# #                         trueDist = TRUE,
# #                         ncex = 1)
# #library(grid)
# #grid.draw(plot)

# ###################################
