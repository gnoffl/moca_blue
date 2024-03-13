#
'DO NOT USE'
'DO NOT USE'
'DO NOT USE'
'SUPER SLOW'
'SUPER SLOW'
FALSE
FALSE
FALSE
USE BCFTOOLS!!!!
  
#Load required packages
library(Biostrings)
library(dplyr)
setwd("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/")
load("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/AthAly.genome.snp.RData")
# Convert the matrix to a data frame
AthAly.genome.snp_df <- as.data.frame(AthAly.genome.snp)
# Create the modifications data frame
modifications <- data.frame(fasta_name = AthAly.genome.snp_df$chr,
                            position_to_be_replaced = AthAly.genome.snp_df$Atpos,
                            character_introduced = AthAly.genome.snp_df$lyr)
# Read the original FASTA file
fasta_file <- readDNAStringSet("./../../ref_seq/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
names_fasta <- names(fasta_file)
#print(names_fasta)
#Translator Table for unusual naming of chromosomes or fasta headers
table_data <- data.frame(
  fasta_name_0 = c("1 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF", 
                   "2 dna:chromosome chromosome:TAIR10:2:1:19698289:1 REF", 
                   "3 dna:chromosome chromosome:TAIR10:3:1:23459830:1 REF", 
                   "4 dna:chromosome chromosome:TAIR10:4:1:18585056:1 REF", 
                   "5 dna:chromosome chromosome:TAIR10:5:1:26975502:1 REF"),
  fasta_name = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
)
############################################
# Function to perform changes on the sequence
######################################################
perform_changes <- function(sequence, input_list) {
  for (i in 1:nrow(input_list)) {
    position <- input_list[i, "position_to_be_replaced"]
    new_char <- input_list[i, "character_introduced"]
    
    # Debugging print statements
    print(paste("Position:", position))
    print(paste("Length of sequence:", nchar(sequence)))
    
    # Check if position is within the length of the sequence
    if (position > nchar(sequence)) {
      stop("Position exceeds length of sequence")
    }
    
    sequence <- paste0(substr(sequence, 1, position - 1), new_char, substr(sequence, position + 1, nchar(sequence)))
  }
  return(sequence)
}
#######################################################
############################################
modifications_0 <- merge(modifications, table_data, by = "fasta_name")
modifications_0 <- modifications_0[, -which(names(modifications_0) == "fasta_name")]
names(modifications_0)[which(names(modifications_0) == "fasta_name_0")] <- "fasta_name"
modifications_0$position_to_be_replaced <- as.numeric(modifications_0$position_to_be_replaced)
############################################
# Assuming modifications_0 is a data frame containing modifications data
# and fasta_file is an XStringSet object containing FASTA sequences
############################################
# Assuming modifications_0 is a data frame containing modifications data
# and fasta_file is an XStringSet object containing FASTA sequences
#Apply modifications to fasta_file
modified_seqs <- lapply(split(modifications_0, modifications_0$fasta_name), function(mod) {
  seq_index <- which(names_fasta == mod$fasta_name[1])
  if (length(seq_index) > 0) {
    perform_changes(as.character(fasta_file[seq_index]), mod)
  } else {
    NA
  }
})

#fasta_file0 <- DNAStringSet(fasta_file)
# Remove non-DNA characters from fasta_file
#filtered_fasta_file <- gsub("[^ACGTacgtNn]", "", fasta_file)
# Convert to DNAStringSet
#fasta_file0 <- DNAStringSet(filtered_fasta_file)
# Write the modified sequences to a new FASTA file
writeXStringSet(fasta_file, "Atha-genom_Alyr-mod_TAIR10.fas")

############################
# Define the input fasta sequence
fasta_sequence <- ">seq1\nATGCCGCTAGCGTGG"

# Define the input list of changes
input_list <- data.frame(sequence = c("seq1", "seq1", "seq1"),
                         position = c(2, 4, 5),
                         new_character = c("C", "A", "A"))
# Function to perform changes on the sequence
######################################################
perform_changes <- function(sequence, input_list) {
  for (i in 1:nrow(input_list)) {
    position <- input_list[i, "position"]
    new_char <- input_list[i, "new_character"]
    sequence <- paste0(substr(sequence, 1, position - 1), new_char, substr(sequence, position + 1, nchar(sequence)))
  }
  return(sequence)
}
#######################################################
# Extract sequence ID and sequence string from the fasta sequence
seq_lines <- unlist(strsplit(fasta_sequence, "\n"))
seq_id <- seq_lines[1]
sequence <- paste(seq_lines[-1], collapse = "")


# Perform changes
new_sequence <- perform_changes(sequence, input_list)

# Output the new fasta sequence with original sequence name
cat(paste0(">", seq_id, "\n", new_sequence, "\n"))