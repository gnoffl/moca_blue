########################## [SETUP] ##############################
working_directory <- "/home/ibg-4/Desktop/Rhome/moca_blue/mo_nom"
#######################
jaspar_file <- "./out/rdf5_epmAsagSW_cwm-motifs.jaspar"
#######################################################
dirpath_out = "../out"
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
########################## [PROCESSING] ##############################
# Set working directory
setwd(working_directory)

# Function to install required packages automatically
install_required_packages <- function() {
  # Check if BiocManager is installed, install if not
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # List of required packages and their sources
  bioc_packages <- c("TFBSTools", "motifStack", "universalmotif", "JASPAR2020")
  
  # Install missing Bioconductor packages
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}

# Install required packages
cat("Checking and installing required packages...\n")
install_required_packages()

# Load required libraries with error handling
tryCatch({
  suppressPackageStartupMessages({
    library(TFBSTools)
    library(motifStack)
    library(universalmotif)
    library(JASPAR2020)
    
    # Use JASPAR2020 database
    jaspar_db <- JASPAR2020
  })
}, error = function(e) {
  cat("Error loading required packages:", e$message, "\n")
  cat("Please ensure the following packages are installed:\n")
  cat("- TFBSTools\n")
  cat("- motifStack\n") 
  cat("- universalmotif\n")
  cat("- JASPAR2020\n")
  stop("Required packages not available")
})

##################################################################################
# Update the file name and path accordingly
output_prefix <- gsub("\\.jaspar$", "", jaspar_file) # For output filenames
##################################################################################

# + + + + + + +#                                                               # + + + + + + + #
output_filename <- file.path(dirpath_out, paste0(basename(output_prefix), "_comparison_JASPAR2020.csv"))
##################################################################################
tryCatch({
  pfm <- read_jaspar(jaspar_file)
  cat("Successfully read", length(pfm), "motifs from JASPAR file\n")
  
  # Check if pfm is valid
  if (length(pfm) == 0) {
    stop("No motifs found in the JASPAR file")
  }
  
  # Print motif classes for debugging
  cat("Motif classes found:", class(pfm), "\n")
  if (is.list(pfm)) {
    cat("First motif class:", class(pfm[[1]]), "\n")
  }
  
}, error = function(e) {
  cat("Error reading JASPAR file:", e$message, "\n")
  stop("Failed to read motif file")
})
pwm_uni0 <- convert_motifs(pfm, class = "TFBSTools-PWMatrix") # Using PWMatrix directly

# Get plant PWMs from the available JASPAR database
jaspar_motifs_plants <- getMatrixSet(jaspar_db,
                                     opts = list(tax_group = 'plants',
                                                 matrixtype = 'PWM'))

###################################################
# Prepare your input motifs for comparison
omni_list <- convert_motifs(pfm, class = "PWMatrix")
####################################################

# Convert JASPAR motifs to PWMatrix for TFBSTools compatibility
jaspar_motifs_plants0 <- convert_motifs(jaspar_motifs_plants, class = "PWMatrix")

###################################################
# Compare each of your motifs against all JASPAR plant motifs and visualize top matches
comparison_results_list <- list()
for (i in seq_along(omni_list)) {
  your_motif <- omni_list[[i]]
  similarities <- compare_motifs(your_motif, jaspar_motifs_plants0,
                                 method = "Pearson", min.overlap = 5) # Adjust as needed
  
  # Sort by similarity score (descending)
  similarities_sorted <- similarities[order(score(similarities), decreasing = TRUE)]
  
  # Extract top matches (e.g., top 3 for visualization and reporting)
  top_n <- 3
  top_matches_indices <- head(which(score(similarities_sorted) > 0), n = top_n)
  
  if (length(top_matches_indices) > 0) {
    top_matches <- similarities_sorted[top_matches_indices]
    
    # Create a data frame for the top matches
    motif_results <- data.frame(
      Your_Motif_ID = ID(your_motif),
      Your_Motif_Name = name(your_motif),
      JASPAR_Motif_ID = ID(top_matches),
      JASPAR_Motif_Name = name(top_matches),
      Similarity_Score = score(top_matches),
      P_value = pval(top_matches),
      E_value = eval(top_matches)
    )
    comparison_results_list[[i]] <- motif_results
    
    # Visualize the comparison using motifStack
    plot_filename <- paste0(output_prefix, "_motif_", ID(your_motif), "_top_matches.pdf")
    pdf(plot_filename, width = 8, height = 6)
    tryCatch({
      motifs_to_plot <- c(list(your_motif), as.list(jaspar_motifs_plants0[names(top_matches)]))
      names(motifs_to_plot) <- c(paste0("YourMotif_", ID(your_motif)), name(top_matches))
      print(motifStack(motifs_to_plot, layout = "stack"))
    }, error = function(e) {
      cat(paste("Error plotting motifs for", ID(your_motif), ":", e$message, "\n"))
    })
    dev.off()
    cat(paste("Motif comparison plot saved to:", plot_filename, "\n"))
    
  } else {
    comparison_results_list[[i]] <- data.frame(
      Your_Motif_ID = ID(your_motif),
      Your_Motif_Name = name(your_motif),
      JASPAR_Motif_ID = NA,
      JASPAR_Motif_Name = NA,
      Similarity_Score = NA,
      P_value = NA,
      E_value = NA
    )
  }
}

# Combine all comparison results into a single data frame
final_comparison_results <- do.call(rbind, comparison_results_list)

###########################################################################
write.table(final_comparison_results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat(paste("Comparison results saved to:", output_filename, "\n"))
###########################################################################