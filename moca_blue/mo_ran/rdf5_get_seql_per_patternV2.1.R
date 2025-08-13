library(rhdf5)
# install.packages("plyr")
library(plyr)
library(tidyr)
###################################################################################
working_directory = "~/Desktop/Rhome/moca_blue/mo_clu"
BASENAME = "SPECIES_MODEL_WHATEVER"
#######################################################
FILE1= "Atha-0h_modisco.hdf5"
###################################################################################
dirpath_out = "../out/ranges"

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 4) {
  stop("Cannot provide more than 4 arguments!\nUsage: <script.R> [working_directory] [BASENAME] [FILE1] [dirpath_out]", call.=FALSE)
}
if (length(args)==4) {
  dirpath_out = args[4]
}
if (length(args)>=3) {
  FILE1 = args[3]
}
if (length(args)>=2) {
  BASENAME = args[2]
}
if (length(args)>=1) {
  working_directory = args[1]
}
# print all arguments
cat("BASENAME:", BASENAME, "\n")
cat("FILE1:", FILE1, "\n")
cat("dirpath_out:", dirpath_out, "\n")
cat("working_directory:", working_directory, "\n")

setwd(working_directory)
# Define the pattern names to iterate over

###################################################
#motifs in metacluster 0 -automatize this step
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#X=20
#motifs in metacluster 1 -automatize this step
#Y=10
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
h5file <- H5Fopen(file.path(FILE1), "H5F_ACC_RDONLY")
h5ls(h5file)
metacluster_group <- h5read(h5file, "metacluster_idx_to_submetacluster_results")
# loop through the metaclusters 0 and 1
for (i in c(0, 1)) {
  metacluster <- metacluster_group[[paste0("metacluster_", i)]]
}
#######################################################
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

##############################################             METACLUSTER_0    ######
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

############################################################################ ######
file_path_out <- file.path(dirpath_out, paste0(BASENAME))
dir.create(dirpath_out, recursive = TRUE, showWarnings = FALSE)
write.csv(df, file = paste0(file_path_out,".txt"), row.names = FALSE)
