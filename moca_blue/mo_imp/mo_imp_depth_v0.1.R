######################
library(rhdf5)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
######################################################
setwd("~/Desktop/Rhome/moca_blue/mo_imp")
###################### Setup for "moca_blue" enviroment
NAME0="rdf5_"
SPEC="Arth"
MODEL="S0" # C0 stand for DeepCistrome version 1 (available at 02-may-2023) Standard conditions
DATE= "20231121"
#######################################################
FILE1= "arabidopsis_scores.h5"
FILE2= "arabidopsis_meta_saliency_info.csv"
#######################################################
dirpath_in1 = "../0MOTIFS/modisco_SSR_nov-edition/"
dirpath_out = "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
file_path1 = file.path(dirpath_in1,FILE1)
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
model_parameter <- read.table(
  file.path(
    dirpath_in1,
    FILE2),
  header=TRUE,
  sep=",")
#######################################################
model_parameter$pred_class <- ifelse(model_parameter$pred_prob >= 0.5, 1, 0)
model_parameter$TRUE_class <- ifelse(model_parameter$pred_class == model_parameter$true_target, TRUE, FALSE)
#######################################################
h5file <- H5Fopen(file_path1, "H5F_ACC_RDONLY")
#h5ls(h5file)
saliency_scores <- h5read(h5file,
                          "contrib_scores")
#######################################################
num_arrays <- dim(saliency_scores)[3]
num_columns <- dim(saliency_scores)[2]
result_matrix <- matrix(0, nrow = num_arrays, ncol = num_columns)
for (i in 1:num_arrays) {
  current_array <- saliency_scores[,,i]
  column_sums <- colSums(current_array)
  result_matrix[i,] <- column_sums
}
imp_scores<-as.data.frame(result_matrix)
##################################################################################
model_parameter$sum_imp_score <- rowSums(imp_scores)
model_parameter$sum_imp_score2 <- ifelse(model_parameter$sum_imp_score >= 0, 1, 0)
model_parameter$TRUE_pred_imp <- ifelse(model_parameter$sum_imp_score2 == model_parameter$true_target, TRUE, FALSE)
gene_imp_scores <- as.data.frame(model_parameter$gene_id)
gene_imp_scores <- cbind(gene_imp_scores, imp_scores)
colnames(gene_imp_scores)[1] <- "loc_ID"
############################################

selected_gene <- imp_scores[1,]
slctd_gene_t0 <- as.data.frame(t(selected_gene))
column_1_sum <- sum(slctd_gene_t0[, 1])
column_1_sum
rooted_sqrd_val <- sqrt(slctd_gene_t0[, 1]^2)
rtd_sqr_sum <- sum(rooted_sqrd_val)
rtd_sqr_sum
slctd_gene_t0$norm_imp <- (slctd_gene_t0[, 1] / column_1_sum)*-1
slctd_gene_t0$norm_prd_imp <- c(sqrt(slctd_gene_t0[, 1]^2)/rtd_sqr_sum)

head(slctd_gene_t0)

plot(1:nrow(slctd_gene_t0), slctd_gene_t0[, 1], type = "l", col = rgb(0, 0, 1, alpha = 0.3), xlab = "Row", ylab = "Values", main = "Line Plot of Columns")
lines(1:nrow(slctd_gene_t0), slctd_gene_t0$norm_imp, col = rgb(1, 0, 0, alpha = 0.3))
lines(1:nrow(slctd_gene_t0), slctd_gene_t0$norm_prd_imp, col = rgb(0, 1, 0, alpha = 0.8))

legend("topright", legend = c("1", "norm_imp", "norm_prd_imp"), col = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5), rgb(0, 1, 0, alpha = 0.5)), lty = 1)




