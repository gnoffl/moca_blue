#'THIS NEEDS AN UPDATE"
# mo_feat-filter.v2.5.R does multiple this. First, it can apply different filters to the parsed mapping data coming from "occ_filter.v1.1.R".
# In the near future, occ_filter [~ 50 lines] and mo_feat-filter [400+  lines] will be joined. In addition, to its filter-fucntion, this script 
# produces different outputs that can be used for downstream analyses. Currently, it features human-readable table files for interpretation, 
# bed-files for visualization and annotation, and output files that are ready to use for functional analyses ( e.g. GO-term enrichment). 
# Accordingly, this script takes over a main role in the moca_blue development. All steps do not require more than standard libraries dplyr 
# and stringr and files that are required for model training or are being produced by upstream procedures of the moca_blue enviroment. 
# 2023-08-08 Dr. Simon M. Zumkeller

#--older documentation --
# /home/ibg-4/Desktop/Rhome/solanum_motifs
# There are four files in total which will'be used to process the data
# There is SolyMSR-TSS_motif_position.csv and SolyMSR-TTS_motif_position.csv
# These two files contain the summary statistics of the motifs and their margins
# Then there is Solanum_lycopersicum.SL3.0.55.chr.gff3 (AND do it for spenn_v2.0_gene_models_annot.gff)
# These files contain information about the genes region freature
# These files are located in the directory spenn_v2.0_gene_models_annot.gff
# Finally, there is the file out_all.txt
# this file is located in the directory /home/ibg-4/Desktop/Rhome/solanum_motifs/out
# It contains which motifs matched what gene region
# All files need to imported into R and combined
# The goal of this script is to characterize genes that have identified motifs in their flanking regions
#gff annotations might be edited eg. sed -i 's/ID=gene:[^A]*AT/ID=gene:AT/g' path/to/file
##############################################################
working_directory = "/home/ibg-4/Desktop/Rhome/moca_blue/mo_proj/"
BASENAME = "SPECIES_MODEL_DATE"
##############################################################
weight_region <- "yes"   # choose between "yes" or "no". If you write yes, an additional filter will be applied to occurences
word_size <- 14  # This is not a real wordsize like in BLAST but it works similar. All matches smaller than this size will be removed.
##############################################################
output_folder <- "../out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
TSS_motifs <- "ZemaS0-TSS_motif_ranges_q1q9.csv"         # These TSS are the summary statistics of the seqlet distribution of the HDF5 file  
TTS_motifs <- "ZemaS0-TTS_motif_ranges_q1q9.csv"         # These TTS are the summary statistics of the seqlet distribution of the HDF5 file
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
annotation_file <- "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gff3" # Genome annotation file  #### CODE HERE EXTRACTS NOT COMPLETE GFFF! NEED TO UPDATE CORRECT !!! ASAP
filtered_mapping_path <- "result_ZemaS0_genZemaB73p0e-4_2024h"    #   These are the filtered results of the the BLAMM output occurence.txt
##############################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 9) {
  stop("Cannot provide more than 9 arguments! Usage: <script> [working_directory] [BASENAME] [weight_region] [word_size] [output_folder] [TSS_motifs] [TTS_motifs] [annotation_file] [filtered_mapping_path]", call.=FALSE)
}
if (length(args)>=9) {
  word_size = as.numeric(args[9])
}
if (length(args)>=8) {
  weight_region = args[8]
}
if (length(args)>=7) {
  filtered_mapping_path = args[7]
}
if (length(args)>=6) {
  annotation_file = args[6]
}
if (length(args)>=5) {
  TTS_motifs = args[5]
}
if (length(args)>=4) {
  TSS_motifs = args[4]
}
if (length(args)>=3) {
  output_folder = args[3]
}
if (length(args)>=2) {
  BASENAME = args[2]
}
if (length(args)>=1) {
  working_directory = args[1]
}
# print all arguments
print(paste("Working directory:", working_directory))
print(paste("Base name:", BASENAME))
print(paste("Weight region:", weight_region))
print(paste("Word size:", word_size))
print(paste("Output folder:", output_folder))
print(paste("TSS file:", TSS_motifs))
print(paste("TTS file:", TTS_motifs))
print(paste("Annotation file:", annotation_file))
print(paste("Filtered mapping results:", filtered_mapping_path))
##############################################################
setwd(working_directory)
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
#   #   #   #   #  #  #   #   #   #  #  #  #  #  #  #  #  #
### FILTERs #### All filters will be applied to the q1q9 output file. 
# For comparison default output shows motif matches in gene flanking regions (1500 kbp up/do window)
# and mima output applies boundaries identified by the pred.-model
# in the up- or downstream region. Acc. to "yes" a motif must appear in more than 10% of cases the up- or downstream region. 
# Available filters for motif_orient: "forward", "reverse", "none"
Filter_motif_orient <- "none"
# Available filters for annot_type : "gene", "CDS", mRNA", "UTR", "none"
Filter_annot_type   <- "gene"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##############################################################
##############################################################
# Load the necessary packages
library(dplyr)
library(stringr)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# # The following initializes usage of Bioc devel
# BiocManager::install(version=3.14, lib="/home/gernot/Rpackages")
# # BiocManager::install("Rhtslib", lib="/home/gernot/Rpackages")
# BiocManager::install("rtracklayer", lib="/home/gernot/Rpackages")
library(rtracklayer)
##############################################################




#THIS MUST BE CHANGED AN ADJUSTED WHEN IT IS USED ON THE OTHER SPECIES !!
# convert_epm <- function(code) {
#   code1 <- substr(code, 1, 12) 
#   code2 <- substr(code, 13, 17)
#   code <- paste0(code1, code2)
#   code <- gsub("Sola", "Soly", code, fixed = TRUE) #ONLY SOLA mispelling
#   code <- gsub("epm", "epm_", code, fixed = TRUE)
#   code <- gsub("__", "_", code, fixed = TRUE)
#   return(code)
# }
convert_epm <- function(code) {
  # Flexible EPM conversion function
  # Handles different species (ARTH, ZEMA, etc.) and dataset types (SSR, MSR, etc.)
  
  # Original formats examples:
  # "epm_SSR_ARTH_S0X0.75dC7K25f_p0m20F_315" -> "epm_Arth_S0X0.75dC7K25f_p0m20"
  # "epm_MSR_ZEMA_S0X0.75dC7K25f_p1m11R_197" -> "epm_Zema_S0X0.75dC7K25f_p1m11"
  
  # Step 1: Remove dataset type prefix (SSR_, MSR_, etc.)
  code <- gsub("epm_(SSR|MSR|[A-Z]+)_", "epm_", code)
  
  # Step 2: Convert species codes to proper case
  # ARTH -> Arth, ZEMA -> Zema, SOLA -> Soly, etc.
  species_mapping <- list(
    "ARTH" = "Arth",
    "ZEMA" = "Zema", 
    "SOLA" = "Soly",
    "SOLY" = "Soly",
    "ATHA" = "Atha"
  )
  
  # Apply species mapping
  for (old_code in names(species_mapping)) {
    pattern <- paste0("_", old_code, "_")
    replacement <- paste0("_", species_mapping[[old_code]], "_")
    code <- gsub(pattern, replacement, code, fixed = TRUE)
  }
  
  # Step 3: Extract pattern base and preserve the actual pattern number
  # Extract everything before the F/R suffix: _p0m20F_315 -> _p0m20
  code <- gsub("_p([0-9]+)m([0-9]+)[FR]_([0-9]+)$", "_p\\1m\\2", code)
  
  return(code)
}






##############################################################
tss_motifs <- read.table(TSS_motifs, header = TRUE
                         , sep=",")
tts_motifs <- read.table(TTS_motifs, header = TRUE
                         , sep=",")

gene_annot <- import(annotation_file)
a<-as.data.frame(gene_annot)
#colnames(a)
#unique(a$type)
# Check if "gene_id" or "ID" is present in column names
target_columns <- c("seqnames",
                    "source",
                    "type",
                    "start",
                    "end",
                    "score",
                    "strand",  
                    "phase")

#print column names of the data frame
if ("gene_id" %in% colnames(a)) {
  target_columns <- c(target_columns, "gene_id")
} else if ("Name" %in% colnames(a)) {
  target_columns <- c(target_columns, "Name")
} else if ("ID" %in% colnames(a)) {
  target_columns <- c(target_columns, "ID")
} else {
  print(colnames(a))
  stop("No suitable identifier column found in the GFF file. Please ensure it contains 'gene_id', 'Name', or 'ID'.")
}

a <- a[, target_columns]
colnames(a) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gene_annot <- a
gene_annotX <- a
############################################ This is superflous and can made shorter ##########
#if (Filter_annot_type == "gene") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("gene", type))
#} else if (Filter_annot_type == "CDS") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("CDS", type))
#} else if (Filter_annot_type == "mRNA") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("mRNA", type))
#} else if (Filter_annot_type == "UTR") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("UTR", type))
#} else if (Filter_annot_type == "five_prime_UTR") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("five_prime_UTR", type))
#} else if (Filter_annot_type == "three_prime_UTR") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("three_prime_UTR", type))
#} else if (Filter_annot_type == "none") {
#  gene_annot <- gene_annot %>%
#    filter(grepl("gene","CDS","mRNA","UTR", type))
#} else {
#  stop("Invalid value for Filter_annot_type")
#}
unique(gene_annot$type)
##############################################################
motif_gene_matches <- read.table(
  file.path(filtered_mapping_path),
  header=TRUE,
  sep="\t")
##############################################################
#Combine DFs by Chr:start-end as loc; Use extracted ranges
#Create loc to idx df




## try different version
gene_annot$loc <- gene_annot$chr
# gene_annot$loc <- paste(
#   gene_annot$chr,
#   ":",
#   gene_annot$start - 1000,
#   "-",
#   gene_annot$end + 1000,
#   sep = "")






#Reduce attributes
#gene_annot$loc_ID <- str_extract(
#  gene_annot$attributes,
#  "(?<=ID=).*?(?=;)")
#gene_annot$loc_ID <- ifelse(
#  grepl(":",
#        gene_annot$loc_ID),
#  sub(".*:", "",
#      gene_annot$loc_ID),
#  gene_annot$loc_ID)
#subset the df
gene_annot$loc_ID <- gene_annot$attributes
gene_annot0 <- gene_annot[, c("loc",
                              "loc_ID",
                              "chr",
                              "type",
                              "start",
                              "end",
                              "strand")]
##################################################
filtered_motif_df <- motif_gene_matches[motif_gene_matches$mstart <= 1500 |
                                          abs(
                                            motif_gene_matches$mstart-(
                                              motif_gene_matches$gene_end - motif_gene_matches$gene_start)+1) <= 1500, ]
##################################################
# Filter the data frame based on Filter_motif_orient variable
if (Filter_motif_orient == "forward") {
  filtered_motif_df <- filtered_motif_df %>%
    filter(grepl("F_", motif))
} else if (Filter_motif_orient == "reverse") {
  filtered_motif_df0 <- filtered_motif_df %>%
    filter(grepl("R_", motif))
} else if (Filter_motif_orient == "none") {
  filtered_motif_df <- filtered_motif_df %>%
    filter(grepl("F_|R_", motif))
} else {
  stop("Invalid value for Filter_motif_orient")
}
##################################################
subset_filmotif_df <- filtered_motif_df[, c("loc",
                                            "motif",
                                            "mstart",
                                            "mend",
                                            "score",
                                            "strand")]
print("sample of subset_filmotif_df")
print(head(subset_filmotif_df, 5))
print("sample of gene_annot0")
print(head(gene_annot0, 5))
merg_df <- merge(gene_annot0,
                 subset_filmotif_df,
                 by = "loc", all.x = TRUE)
merg_df <- na.omit(merg_df)
# Internal check for the number of rows
if (nrow(merg_df) < 2) {
  error_message <- "Please check identifiers chr#:start-end for fasta and gff input"
  stop(error_message)
}
################################################## RC RC RC RC RC CALCULATION OF GENOMIC AND EXTRACTED POSITIONS
merg_df$gen_start = merg_df$start - 1000
merg_df$gen_end = merg_df$end + 1000
merg_df$gen_mstart = merg_df$mstart + merg_df$gen_start # !!!!
#merg_df$gen_mstart <- ifelse(merg_df$strand.y =="-", merg_df$gen_mstart+1, merg_df$gen_mstart)
merg_df$gen_mend = merg_df$mend + merg_df$gen_start
#merg_df$gen_mend <- ifelse(merg_df$strand.y =="-", merg_df$gen_mend+1, merg_df$gen_mend)
merg_dfA <- merg_df
merg_dfA$dis_gms_gst = abs(merg_dfA$gen_mstart - merg_dfA$gen_start)
merg_dfA$dis_gms_ge = abs(merg_dfA$gen_mstart - merg_dfA$gen_end)
merg_dfA$dis_gme_gst = abs(merg_dfA$gen_mend - merg_dfA$gen_start)
merg_dfA$dis_gme_ge = abs(merg_dfA$gen_mend - merg_dfA$gen_end)
merg_dfA$dist_gen_mstart <- pmin(merg_dfA$dis_gms_gst, merg_dfA$dis_gms_ge)
merg_dfA$dist_gen_mend <- pmin(merg_dfA$dis_gme_gst, merg_dfA$dis_gme_ge)
merg_dfA$dist_gen_mo_len <- 1+abs(merg_dfA$dist_gen_mstart - merg_dfA$dist_gen_mend)
merg_dfB <- merg_dfA[merg_dfA$dist_gen_mo_len == word_size, ]
merg_dfB <- merg_dfB %>%
  select(chr, loc, loc_ID, type, start, end, strand.x, motif, dist_gen_mstart, dist_gen_mend, gen_mstart, gen_mend, strand.y)
merg_dfB <-merg_dfB %>%
  mutate(
    dist_gen_mstart = ifelse(strand.y == "-", dist_gen_mstart + 1, dist_gen_mstart),
    dist_gen_mend = ifelse(strand.y == "-", dist_gen_mend + 1, dist_gen_mend),
    gen_mstart = ifelse(strand.y == "-", gen_mstart + 1, gen_mstart),
    gen_mend = ifelse(strand.y == "-", gen_mend + 1, gen_mend)
  )
merg_dfB$region <- ifelse(
  merg_dfB$dist_gen_mstart  < 
    merg_dfB$dist_gen_mend,
  "upstream",
  "downstream")
############################################## IN PROGRESSS ###
################################################################
#merg_df$gen_mstart <- merg_df$mstart -1000 + merg_df$start
#merg_df$dist_gen_start <- abs(merg_df$gen_mstart - merg_df$start +1000) 
#merg_df$dist_gen_end <- abs(merg_df$end + 1000 - merg_df$gen_mstart )
#############################################################################
#merg_df<- subset(merg_df, score >= "20")
############################################################################# THIS STEP IS NOT TRIVIAL
########################### TRANSCRIPTS THAT ARE SMALLER THAN 1K CAN HAVE MATCHES IN UP AND DOWN STREAM
############# REGIONS. THERE IS JUST ONE MATCH BUT IT IS EQUALLY REPRESENTED IN BOTH DATA SETS.
#### FILTER BY MOTIF PREFERENCES ########## MAY INCLUDE FLAG THAT MATCH APPEARS IN UP AND DOWN - LATER 
# Subsetting columns from merg_df to create merg_df1 without dist_gen_start ###########################
merg_df1 <- merg_dfB %>%
  select(chr, loc, loc_ID, type, start, end, strand.x, motif, dist_gen_mstart, gen_mstart, gen_mend, strand.y, region)
colnames(merg_df1)[colnames(merg_df1) == "dist_gen_mstart"] <- "dist_transc_border"
merg_df00 <- distinct(merg_df1)
merg_df00 <- merg_df00[merg_df00$dist_transc_border <= 1500, ]
#merg_df1 <- merg_dfB[, !(colnames(merg_dfB) %in% c("dist_gen_mstart"))]
#merg_df2 <- merg_dfB[, !(colnames(merg_dfB) %in% c("dist_gen_mend"))]
#colnames(merg_df1)[colnames(merg_df1) == "dist_gen_mstart"] <- "dist_transc_border"
#colnames(merg_df2)[colnames(merg_df2) == "dist_gen_mstart"] <- "dist_transc_border"
#colnames(merg_df2)[colnames(merg_df2) == "dist_gen_mend"] <- "dist_transc_border"
#merg_df01 <- rbind(merg_df1, merg_df2)
#############################################################################
# Invert region based on strand.x column
merg_df00$region <- ifelse(
  merg_df00$strand.x == "-",
  ifelse(
    merg_df00$region == "upstream",
    "downstream",
    "upstream"),
  merg_df00$region)
##### CREATE ONE OUTPUT FILE THAT CONTAINS UNSPECIFIC MATCHES ###############
#############################################################################
#############################################################################
##### START TO FILTER USING THE MOTIF SEQLET SUMMARY STATISTICS #############
merg_df00$epm <- convert_epm(merg_df00$motif)
## #### ####### #################### ######################## NOW COMPLICATED STUFF
up_merg_df00 <- merg_df00[merg_df00$region == "upstream", ]
#subset_up_merg <- up_merg_df00[, c(1, 14, 15, 16)]
subset_tss_m <- tss_motifs[, c(1, 2, 3, 7, 8)]
merg_up_mo <- merge(up_merg_df00,
                    subset_tss_m,
                    by = "epm", all.x = TRUE)
merg_up_mo0 <- na.omit(merg_up_mo)
# Internal check for the number of rows
if (nrow(merg_up_mo0) < 2) {
  error_message <- "Please check identifiers in epm for TSS/TTS ranges.csv and occurences.txt"
  stop(error_message)
}
merg_up_mima_mo0 <- merg_up_mo0 %>%
  filter(dist_transc_border >= min & dist_transc_border <= max)
merg_up_q10_q90_mo0 <- merg_up_mo0 %>% 
  filter(dist_transc_border >= q10 & dist_transc_border <= q90)
do_merg_df00 <- merg_df00[merg_df00$region == "downstream", ]
subset_tts_m <- tts_motifs[, c(1, 2, 3, 7, 8)]
#subset_tts_m$min <- subset_tts_m$min -1520
#subset_tts_m$max <- subset_tts_m$max -1520
#subset_tts_m$q10 <- subset_tts_m$q10 -1520
#subset_tts_m$q90 <- subset_tts_m$q90 -1520
merg_do_mo <- merge(do_merg_df00,
                    subset_tts_m,
                    by = "epm", all.x = TRUE)
merg_do_mo0 <- na.omit(merg_do_mo) #! ???????????????????????????????????
# print("Sample downstream regions:")
# print(head(do_merg_df00, 5))
# print("Sample merged downstream motif data before removing NAs:")
# print(head(merg_do_mo, 5))
# print("Sample merged downstream motif data after removing NAs:")
# print(head(merg_do_mo0, 5))
# Internal check for the number of rows
if (nrow(merg_do_mo0) < 2) {
  error_message <- "Please check identifiers in epm for TSS/TTS ranges.csv and occurences.txt"
  stop(error_message)
}

### script crashes here:

merg_do_mima_mo0 <- merg_do_mo0 %>%
  filter(dist_transc_border >= min & dist_transc_border <= max)
merg_do_q10_q90_mo0 <- merg_do_mo0 %>% 
  filter(dist_transc_border >= q10 & dist_transc_border <= q90)

a_mima_df01 <- rbind(merg_up_mima_mo0,
                     merg_do_mima_mo0)
a_q10q90_df01 <- rbind(merg_up_q10_q90_mo0,
                       merg_do_q10_q90_mo0)

a0_DF_df <- merg_df00
a0_mima_df <- a_mima_df01
a0_q1q9_df <- a_q10q90_df01

#a0_DF_df <- merg_df00[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
#a0_mima_df <- a_mima_df01[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
#a0_q1q9_df <- a_q10q90_df01[, c(2, 3, 5, 6, 7,8,15,16,9,14,13,1,11,12)]
############################################################################## APPLY WEIGHT FILTER
if (weight_region == "yes") {
  tss_motifs$weighted_region <- NA
  tts_motifs$weighted_region <- NA
  for (i in tss_motifs$epm) {
    tts_row <- tts_motifs[tts_motifs$epm == i, ]
    # Check if matching TTS row exists
    if (nrow(tts_row) == 0) {
      # No matching TTS motif found, skip this iteration
      cat("Warning: No TTS match found for EPM:", i, "\n")
      next
    }

    c_sum <- sum(tss_motifs[tss_motifs$epm == i, "number"], tts_row$number)
    threshold <- 0.2 * c_sum                                             ################ change weight here here here here
    if (tss_motifs[tss_motifs$epm == i, "number"] < threshold) {
      tss_motifs[tss_motifs$epm == i, "weighted_region"] <- 1
    } else {
      tss_motifs[tss_motifs$epm == i, "weighted_region"] <- 0
    }
    if (tts_row$number < threshold) {
      tts_motifs[tts_motifs$epm == i, "weighted_region"] <- 1
    } else {
      tts_motifs[tts_motifs$epm == i, "weighted_region"] <- 0
    }
  }
  tss_motifs_ss <-  tss_motifs[,c(1,14)]
  tss_motifs_ss$region <-  c("upstream")
  tts_motifs_ss <-  tts_motifs[,c(1,14)]
  tts_motifs_ss$region <-  c("downstream")
  region_weights <- rbind(tss_motifs_ss, tts_motifs_ss)
  region_weights <- region_weights %>%
    filter(weighted_region != 0)
  region_weights0 <- region_weights[,c(1,3)]
  a0_q1q9_df <- a0_q1q9_df %>%
    anti_join(region_weights0, by = c("epm", "region"))
} else {
  a0_q1q9_df <- a_q10q90_df01[, c(2, 3, 5, 6, 7, 8, 15, 16, 9, 14, 13, 1, 11, 12)]
}
################################################## (similar to mo_finder line 139, expect for feature filter)
#Create loc to idx df
#gene_annotX$loc <- paste(
#  gene_annotX$chr,
#  ":",
#  gene_annotX$start - 1000,
#  "-",
#  gene_annotX$end + 1000,
#  sep = "")
#Reduce attributes
#gene_annotX$loc_ID <- str_extract(
#  gene_annotX$attributes,
#  "(?<=ID=).*?(?=;)")
#gene_annotX$loc_ID <- ifelse(
#  grepl(":",
#        gene_annotX$loc_ID),
#  sub(".*:", "",
#      gene_annotX$loc_ID),
#  gene_annotX$loc_ID)
#subset the df
#colnames(gene_annotX)
#colnames(gene_annotX) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "loc_ID", "loc" )
#gene_annotX <- gene_annotX[, c("loc",
#                               "loc_ID",
#                               "chr",
#                               "type",
 #                              "start",
#                               "end",
#                               "strand")]
#x <- gene_annotX[gene_annotX$type != "gene", ]
#unique(x$type)
#x <- x[, c("loc_ID",
#           "loc",
#           "type",
#           "start",
#           "end")]
#x$loc_ID <- sub("\\..*", "", x$loc_ID)
#y <- gene_annotX[gene_annotX$type == "gene", ]
############################################################################ 'HIER WURM WFÜR MEHR ALS GENERISCHE FEature
#CREATE THE FEATURE "untranscribed"
#gene_annotU <- gene_annot0
#gene_annotU$type.a <- c("untranscr")
#gene_annotU$start.a <- gene_annotU$start -1000
#gene_annotU$end.a <- gene_annotU$start -1
#gene_annotU$start.b <- gene_annotU$end +1
#gene_annotU$end.b <- gene_annotU$end +1000
#gene_annotU1 <- gene_annotU[, c("loc","loc_ID","chr","type.a","start.a","end.a","strand")]
#colnames(gene_annotU1) <- c("loc","loc_ID","chr","type","start","end","strand")
#gene_annotU2 <- gene_annotU[, c("loc","loc_ID","chr","type.a","start.b","end.b","strand")]
#colnames(gene_annotU2) <- c("loc","loc_ID","chr","type","start","end","strand")
#gene_annotU12 <- rbind(gene_annotU1, gene_annotU2)
#z <- gene_annotU12[, c("loc_ID","loc","type","start","end")]
#x0 <- rbind(x,z)
#x0$loc_ID <- gsub("[-']+", "", x0$loc_ID)
#x0$loc_ID <- gsub("\\..*", "", x0$loc_ID)
#y$loc_ID <- gsub("[-']+", "", y$loc_ID)
#y$loc_ID <- gsub("\\..*", "", y$loc_ID)
#unique(x0$type)
###########################################################################
#feat_annot0<- merge(y, x0, by = "loc")
#unique(feat_annot0$type.x)
#feat_annot1<-feat_annot0 # start change here
#feat_annot2 <- feat_annot1[, c("loc_ID.x","loc","type.y","start.y","end.y")]
#colnames(feat_annot2) <- c("loc_ID","loc","type.y","start.y","end.y")
#feat_annot2 <- na.omit(feat_annot2)
#feat_annot2 <- feat_annot2[!duplicated(feat_annot2),]
#a_mima_df01 <- na.omit(a_mima_df01)
#a_mima_df01 <- a_mima_df01[!duplicated(a_mima_df01),] #### CHECK !!!! SCRIPT RUNS TO THIS LINE (2023-08-02)
##########################################################
#a_mima_df01$loc_ID <- gsub("[-']+", "", a_mima_df01$loc_ID)
#a_mima_df01$loc_ID <- gsub("\\..*", "", a_mima_df01$loc_ID)
#feat_anno_mima <- merge (a_mima_df01,
#                         feat_annot2,by = "loc")
#feat_anno_mima0 <- na.omit(feat_anno_mima)           #### CHECK !!!! SCRIPT RUNS TO THIS LINE (2023-08-02)
#feat_anno_mima01 <- subset(feat_anno_mima0, gen_mstart > start.y & gen_mstart < end.y)
#feat_anno_mima01 <- feat_anno_mima01[!duplicated(feat_anno_mima01),] #RUNS!
#unique(a_mima_df01$type)
feat_anno_mima<- a0_q1q9_df
feat_anno_mima$type.y <- ifelse(feat_anno_mima$dist_transc_border >= 1000, "transcribed", "non-transcribed")
feat_anno_mima01 <- feat_anno_mima[, c("loc",
                                         "loc_ID",
                                         "type",
                                         "start",
                                         "end",
                                         "strand.x",
                                         "motif",
                                         "gen_mstart",
                                         "dist_transc_border",
                                         "strand.y",
                                         "region",
                                         "type.y",
                                         "min",
                                         "max",
                                         "q10",
                                         "q90")]
feat_anno_mima01 <- feat_anno_mima01[!duplicated(feat_anno_mima01),]
feat_anno_mima01$epm <- convert_epm(feat_anno_mima01$motif)
#feat_anno_mimaW01 <- feat_anno_mima01[!(feat_anno_mima01$epm %in% region_weights0$epm & feat_anno_mima01$region %in% region_weights0$region), ]
############################################################## CONTINUE HERE --- CREATE BED FILES

a0_mima_df$chr <- sub(":.*", "", a0_mima_df$loc)
a0_q1q9_df$chr <- sub(":.*", "", a0_q1q9_df$loc)
a0_DF_bed <- a0_DF_df[, c("chr","gen_mstart","gen_mend","motif","strand.y")]
a0_mima_bed <- a0_mima_df[, c("chr","gen_mstart","gen_mend","motif","strand.y")]
a0_q1q9_bed <- a0_q1q9_df[, c("chr","gen_mstart","gen_mend","motif","strand.y")]
colnames(a0_DF_bed) <- c("chrom","chromStart","chromEnd","name","strand")
colnames(a0_mima_bed) <- c("chrom","chromStart","chromEnd","name","strand")
colnames(a0_q1q9_bed) <- c("chrom","chromStart","chromEnd","name","strand")
a0_DF_bed$thickStart <- c(a0_DF_bed$chromStart)
a0_DF_bed$thickStart <- c(a0_DF_bed$chromEnd)
a0_DF_bed$itemRgb <- c("255,0,127")
a0_mima_bed$thickStart <- c(a0_mima_bed$chromStart)
a0_mima_bed$thickStart <- c(a0_mima_bed$chromEnd)
a0_mima_bed$itemRgb <- c("255,0,127")
a0_q1q9_bed$thickStart <- c(a0_q1q9_bed$chromStart)
a0_q1q9_bed$thickStart <- c(a0_q1q9_bed$chromEnd)
a0_q1q9_bed$itemRgb <- c("255,0,127")

a0_q1q9_bed <- unique(a0_q1q9_bed)
a0_mima_bed <- unique(a0_mima_bed)
############################################################## GENERATE OUTPUT
file_path_out <- file.path(output_folder, paste0(BASENAME, "_",Filter_annot_type,"_", Filter_motif_orient))
# make folders for output files
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
#write.csv(a0_DF_df, file=paste0(file_path_out,
#                                "-def.csv"), row.names=FALSE)
write.csv(a0_mima_df, file=paste0(file_path_out,
                                  "-mima.csv"), row.names=FALSE)
write.csv(a0_q1q9_df, file=paste0(file_path_out,
                                  "-q1q9.csv"), row.names=FALSE)
#write.table(a0_DF_bed, file=paste0(file_path_out,"-def.bed"), row.names=FALSE, sep = "\t", quote = FALSE)
write.table(a0_mima_bed, file = paste0(file_path_out, "-mima.bed"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(a0_q1q9_bed, file=paste0(file_path_out, "-q1q9.bed"), row.names=FALSE, sep = "\t", quote = FALSE)
write.table(feat_anno_mima01, file=paste0(file_path_out,"feat_mima_weight.csv"), row.names = FALSE, sep = "\t",)
print(paste0(file_path_out,"-def.csv"))
print(paste0(file_path_out,"-mima.csv"))
print(paste0(file_path_out,"-q1q9.csv"))
##############################################################
