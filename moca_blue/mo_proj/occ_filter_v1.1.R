#This script filters the occurences of motifs in a genome (mapped with BLAMM)
#Motifs must lie in a range of 1500 bp of gene start and end, respectively
#Before the occurences can be filtered regarding the motif preferences the genes orientations must be determined
# ATTENTION! The occurence file from BLAMM can be very large.
# To keep computational power low it is highly recommended to split the occurence file into smaller ones (with 1 billion lines)
#############################################################
######################################### USE THIS AS BASH ##
#project="Soly_blamm_20230417"
#inputFile="occurrences_0.00001.txt"
#outputSize="1000000"
#split -l $outputSize --numeric-suffixes $inputFile smallfile
#mkdir occ$project
#mv smallfile* occ$project/
############################################################# ITS WORKING !
working_directory = "."
file_paths = "moca_blue/out/ara_msr_max_corrected/mapped_motifs/occurrences.txt"
output_file <- "moca_blue/out/ara_msr_max_corrected/filtered/ara_msr_max_corrected_filtered.tsv"
dirpath <- "./occSolyM014sols_e3-cwm_20231110/"
# file_paths  <- list.files(dirpath, pattern = "^smallfile*", full.names = TRUE)
#############################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 4) {
  stop("Cannot provide more than 4 arguments! Usage: <script> [<working_directory>] [<file_paths>] [<output_file>] [<dirpath>]", call.=FALSE)
}
if (length(args)>=4) {
  dirpath = args[4]
}
if (length(args)>=3) {
  output_file = args[3]
}
if (length(args)>=2) {
  file_paths = args[2]
}
if (length(args)>=1) {
  working_directory = args[1]
}

if (file.exists(file_paths)){
  file_paths = c(file_paths)
} else if (dir.exists(dirpath)) {
  # list all files in the directory that match the pattern
  file_paths <- list.files(dirpath, pattern = "^smallfile.*\\.txt$", full.names = TRUE)
} else {
  stop("No valid file or directory found.")
}

# print all arguments
print(paste("directory path:", dirpath))
print(paste("Output file:", output_file))
print(paste("working directory:", working_directory))
print(paste("File paths:", paste(file_paths, collapse = ", ")))
############################################################# 
library(tidyr)
library(dplyr)
library(readr)
library(magrittr)
############################################################# 
setwd(working_directory)
# initialize an empty data frame to store the results
combined_df <- data.frame()
# loop over all files and append the results to the combined data frame
for (filepath in file_paths) {
  # read the file and skip lines with less than 9 elements
  occurrence_df <- read.table(filepath, header = FALSE, fill = TRUE)
  # select rows with the specified pattern and update the column names
#  occurrence_df <- subset(occurrence_df, grepl("^epmSola", V3))
  colnames(occurrence_df) <- c("loc",
                               "source",
                               "motif",
                               "mstart",
                               "mend",
                               "score", 
                               "strand",
                               "V7",
                               "V8")
  









  # split the "loc" column into "chr", "gene_start", and "gene_end"
 # selected_df <- separate(occurrence_df, loc, into = c("chr", "gene_start", "gene_end"), sep = "[:-]", fill = "left")
  # selected_df <- separate(occurrence_df,
  #                         loc,
  #                         into = c("chr", "gene_loc"),
  #                         sep = ":",
  #                         fill = "left")  %>%
  #   separate(
  #     gene_loc, into = c("gene_start", "gene_end"),
  #     sep = "-",
  #     fill = "left")
  
  # # combine "chr", "gene_start", and "gene_end" into a new "loc" column
  # selected_df$loc <- paste(selected_df$chr,
  #                          selected_df$gene_start,
  #                          sep = ":")
  # selected_df$loc <- ifelse(
  #   is.na(
  #     selected_df$gene_end),
  #   selected_df$loc,
  #   paste(
  #     selected_df$loc,
  #     selected_df$gene_end,
  #     sep = "-"))
selected_df <- separate(occurrence_df,
                          loc,
                          into = c("chr", "gene_loc"),
                          sep = ":",
                          fill = "left")  %>%
    separate(
      gene_loc, into = c("gene_start", "gene_end"),
      sep = "-",
      fill = "left")
  
# combine "chr", "gene_start", and "gene_end" into a new "loc" column
selected_df$loc <- ifelse(
  any(is.na(selected_df$chr), is.na(selected_df$gene_start)),
  selected_df$gene_end,
  paste(selected_df$chr, selected_df$gene_start, sep = ":")
)

selected_df$loc <- ifelse(
  is.na(selected_df$gene_start),
  selected_df$loc,
  paste(
    selected_df$loc,
    selected_df$gene_end,
    sep = "-"))
selected_df$gene_start = 1001
selected_df$gene_end = 2020













  
  # filter rows based on the specified condition
  filtered_df <- selected_df %>%
    mutate(across(c(mstart, gene_start, gene_end),
                  as.numeric)) %>% # convert columns to numeric
    filter(mstart <= 1500 | abs(mstart - (gene_end-gene_start +1)) <= 1500)
  
  # append the filtered data frame to the combined data frame
  combined_df <- rbind(combined_df, filtered_df)
}
# create folders to save the output
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
# write the combined data frame to a file
write.table(combined_df, 
            output_file, 
            sep = "\t",
            row.names = FALSE)
