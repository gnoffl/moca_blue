if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.13")
BiocManager::install("topGO")



##############################################################
PROJECT <- "SolyMSR_Sope_test"
SPEC <- "Spenn"
MODEL <- "MSR"
DATE <- "20230530"
DATA_ORIGIN <- "motif_matches"
##############################################################
##############################################################
dirpath_1 <- "../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "SolyMSR_on_Spe-ch01-0e3_gene_none20230530feat_mima_q1q9.csv"
file3 <- "mapman_sopen.txt"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path_out <- file.path(dirpath_2, paste0(DATE,"_",PROJECT,"_mo-feat"))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
FILTER<- "q1q9"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep=",")

mapman <- read.table(file3,
                     header=TRUE,
                     sep="\t", quote = "")

########### ############## ############ #####################

colnames(mapman)[colnames(mapman)=="IDENTIFIER"]<-"loc_ID"

mm0$loc_ID <- tolower(mm0$loc_ID)

mapman$loc_ID <- tolower(mapman$loc_ID)
mapman$loc_ID <- gsub("[-']+", "", mapman$loc_ID)
mapman$loc_ID <- gsub("\\..*", "", mapman$loc_ID) #CAREFULL WITH THE DOTS

mapman0 <- mapman[, c("loc_ID",
                      "BINCODE",
                      "NAME")]

mm0_mapman <- merge(mm0, mapman0, by= "loc_ID")

########### ############## ############ ##########

cont_table_A <- table(mm0_mapman$motif, mm0_mapman$BINCODE)

testA =as.data.frame(cont_table_A)
tfestA = testA %>%
  pivot_wider(names_from = Var2, values_from = Freq)

# # # # # # # # # # #  # # # # # # 

cont_table_B <- table(mm0_mapman$motif, mm0_mapman$type.y)

testB =as.data.frame(cont_table_B)
tfestB = testB %>%
  pivot_wider(names_from = Var2, values_from = Freq)

##########################################################################################
##########################################################################################
##########################################################################################

# Check if "mRNA" and "untranscr" columns are present in the data frame
if ("mRNA" %in% colnames(tfestB) && "untranscr" %in% colnames(tfestB)) {
  tfestB$transcr_class <- apply(tfestB[, c("mRNA", "untranscr")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transcr_pref <- ifelse(tfestB$mRNA < tfestB$untranscr, "untranscribed", "transcribed")
} else {
  print("Error: 'mRNA' and/or 'untranscr' columns are not present in the data frame.")
}
##########################################################################################
# Check if "exon" and "intron" columns are present in the data frame
if ("exon" %in% colnames(tfestB) && "intron" %in% colnames(tfestB)) {
  tfestB$feature_class <- apply(tfestB[, c("exon", "intron")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$feature_pref <- ifelse(tfestB$exon < tfestB$intron, "intronic", "exonic")
} else {
  print("Error: 'exon' and/or 'intron' columns are not present in the data frame.")
}
########################################################################################## Comparison might be unfair UTR/CDS better
# Check if "CDS" and "intron" columns are present in the data frame
if ("CDS" %in% colnames(tfestB) && "intron" %in% colnames(tfestB)) {
  tfestB$transl_class <- apply(tfestB[, c("CDS", "intron")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transl_pref <- ifelse(tfestB$CDS < tfestB$intron, "intronic", "codogenic")
} else {
  # Handle the case when "CDS" and/or "intron" columns are not present
  # Print an error message or perform alternative actions
  print("Error: 'CDS' and/or 'intron' columns are not present in the data frame.")
}
##########################################################################################
# Check if "CDS" and "intron" columns are present in the data frame
if ("CDS" %in% colnames(tfestB) && "UTR" %in% colnames(tfestB)) {
  tfestB$transl_class <- apply(tfestB[, c("CDS", "UTR")], 1, function(row) {
    result <- chisq.test(row, p = c(0.5, 0.5))
    result$p.value
  })
  tfestB$transl_pref <- ifelse(tfestB$CDS < tfestB$intron, "UTR", "codogenic")
} else {
  # Handle the case when "CDS" and/or "intron" columns are not present
  # Print an error message or perform alternative actions
  print("Error: 'CDS' and/or 'UTR' columns are not present in the data frame.")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##########################################################################################
##########################################################################################

write.csv(tfestB, file=paste0(file_path_out,
                                "-features.csv"), row.names=FALSE)
##########################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Remove columns with sum less than 5
########### ############## ############ ##########

head(tfestB)
