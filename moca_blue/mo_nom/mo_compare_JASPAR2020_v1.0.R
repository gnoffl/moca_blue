########################## [SETUP] ##############################
working_directory <- "/home/ibg-4/Desktop/Rhome/moca_blue/mo_nom"
#######################
jaspar_file <- "./out/rdf5_ArthD0D6_cwm-motifs.jaspar"
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
# Load required libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("JASPAR2020")
library(TFBSTools)
library(motifStack)
library(universalmotif)
library(JASPAR2020)
##################################################################################
# + + + + + + +#                                                 # + + + + + + + #
# split the extension of the file name and add suffix "_comparison_JASPAR2020.csv"
file_name_without_ext <- tools::file_path_sans_ext(basename(jaspar_file))
output_filename <- file.path(dirpath_out, paste0(file_name_without_ext, "_comparison_JASPAR2020.csv"))
##################################################################################
pfm <- read_jaspar(jaspar_file)
pwm_uni0<-convert_motifs(pfm
  , class = "TFBSTools-PWMatrix")
jaspar_motifs_plants <- getMatrixSet(JASPAR2020,
                                     opts = list(tax_group='plants',
                                                 matrixtype='PWM'))
###################################################
similarities_list <- vector("list", length = length(pwm_uni0))
####################################################
jaspar_motifs_plants0 <-convert_motifs(jaspar_motifs_plants
                                       , class = "TFBSTools-PWMatrix")
###################################################
omni_list<- c(pwm_uni0, jaspar_motifs_plants0)
#omni_comp<- compare_motifs(omni_list, method = "PCC")
#semi_comp <- omni_comp[, !grepl("^epm", colnames(omni_comp))]
#rows_to_keep <- grepl("^epm", rownames(omni_comp))
#semi_comp0f <- semi_comp[rows_to_keep, ]
####################################################
omni_pvlist <- list()
for (i in 1:length(pwm_uni0)) {
  omni_pvlist[[i]] <- compare_motifs(omni_list, i, method = "SW")
}
####  #################################     #########
omni_comparsion_results <- data.frame()
for (i in 1:length(omni_pvlist)) {
  
subject <- omni_pvlist[[i]]@listData[["subject"]]
target <- omni_pvlist[[i]]@listData[["target"]]
score <- omni_pvlist[[i]]@listData[["score"]]
Pval <- omni_pvlist[[i]]@listData[["Pval"]]
Eval <- omni_pvlist[[i]]@listData[["Eval"]]

epm_positions <- which(grepl("^epm", target))

target0 <- target[-(epm_positions)]
score0 <- score[-(epm_positions)]
Pval0 <- Pval[-(epm_positions)]
Eval0 <- Eval[-(epm_positions)]

get <- which.max(score0)

subj <- subject[get]
targ <- target0[get]
scor <- score0[get]
pval <- Pval0[get]
eval <- Eval0[get]

result <- as.data.frame(cbind(subj, targ, scor, pval, eval))
omni_comparsion_results <- rbind(omni_comparsion_results, result)
}
# ## # ## # ## # ## # ## #
ID_df <- data.frame()
for (i in 1:length(jaspar_motifs_plants0)) {
   ID <-  jaspar_motifs_plants0[[i]]@ID
   targ <-  jaspar_motifs_plants0[[i]]@name
   
   result0 <- as.data.frame(cbind(ID, targ))
   ID_df <- rbind(ID_df, result0)
}
###########################################################################
omni_comparsion_results0 <- merge(ID_df, omni_comparsion_results, by = "targ")
omni_comparsion_results0 <- omni_comparsion_results0[, c(3,2,1,4,5,6)]

write.table(omni_comparsion_results0, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
###########################################################################