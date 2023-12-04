library(caTools)
library(caret)
library(randomForest)
library(dplyr)
library(tidyr)
library(gplots)
##############################################################
PROJECT <- "Arth-0e3-cwm-W2q1q9-2"
SPEC <- "Arth"
MODEL <- "S0"
DATE <- "20230823"
##############################################################
dirpath_1 <- "../../ref_seq"
dirpath_2 <- "./out"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
file1 <- "20230823ArthS0Arth-0e3-cwm-W2_gene_nonefeat_mima_weight.csv" #file must contain "feat_"
file2 <- "Arth_S0_predictions.csv" # file with predicted probabilies for expression 0-1
file4 <- "Arabidopsis_thaliana_TPMs-peleke-etal2023.csv"   # file with measurement of models (TPM, quartile classes 0,1,2)
file_path_in_file4 <- file.path(dirpath_1, file4)
file_path_out <- file.path(dirpath_2, paste0(DATE,"_",PROJECT,"_mo-feat"))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
FILTER<- "q1q9" # When q1q9 filter should not be applied, write "NONE", e.g. for min.max filter
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#head(mm0)###############################################################################################
convert_epm <- function(code) {
  code1 <- substr(code, 1, 12) 
  code2 <- substr(code, 13, 18)
  code <- paste0(code1, code2)
  code <- gsub("Sola", "Soly", code, fixed = TRUE) #ONLY SOLA mispelling
  code <- gsub("epm", "epm_", code, fixed = TRUE)
  code <- gsub("__", "_", code, fixed = TRUE)
  return(code)
}
#######################################################################################################
mm0 <- read.table(
  file.path(
    dirpath_2,
    file1),
  header=TRUE,
  sep="\t")
########### ############## ############ #####################
colnames(mm0)

if (FILTER != "NONE") {
  mm0 <- mm0 %>%
    filter(dist_transc_border >= q10 & dist_transc_border <= q90)
} else {
  mm0 <- mm0 %>%
    filter(dist_transc_border >= min & dist_transc_border <= max)
}
##############################################################
head(mm0)
mm0 <- mm0[, c("loc", "loc_ID","type","start","end","strand.x","motif","gen_mstart","dist_transc_border","strand.y","region","type.y")]
colnames(mm0) <- c("loc","loc_ID","type","start","end","loc_strand","motif","gen_mstart","dist_transc_border","mot_strand","region","type.y")
model1 <- read.table(
  file_path_in_file4,
  header=TRUE,
  sep=",")
if (!file.exists(file2)) {
  mm0 <- read.table(file.path(dirpath_2, file1), header=TRUE, sep=",")
  unique_loc_ID <- unique(mm0$loc_ID)  # Remove duplicates
  model0 <- data.frame(loc_ID = unique_loc_ID, prob = sample(c(0, 1), size = length(unique_loc_ID), replace = TRUE))
  print("file2 missing")
  file2_state <- c("FALSE")
} else {
  model0 <- read.table(file2, header=TRUE, sep=",")
  colnames(model0) <- c("loc_ID", "prob")
  print("file2 exists")
  file2_state <- c("TRUE")
}
model0 <- model0[, c(1, 2)]
colnames(model0) <- c("loc_ID", "pred_class")
merged_df <- merge(mm0, model0, 
                   by = "loc_ID",
                   all.x = TRUE,
                   all.y = TRUE)
merged_df <- na.omit(merged_df)
ss_model2 <- model1[, c("gene_id", "true_target")]
colnames(ss_model2) <- c("loc_ID", "expr_class")
ss_model2 <- subset(ss_model2, expr_class !=  2)   
ss_model2$loc_ID <- tolower(ss_model2$loc_ID)
merged_df$loc_ID <- tolower(merged_df$loc_ID)
merg_df01mm<- merge(
  merged_df, ss_model2,
  by = c("loc_ID"),
  all = FALSE, ignore.case = TRUE)
merg_df01mm <- na.omit(merg_df01mm)
merg_df01mm <- merg_df01mm[!duplicated(merg_df01mm),]
unique(merg_df01mm)
merg_df01mm$epm <- convert_epm(merg_df01mm$motif)
testset <- merg_df01mm
testset2 <- testset
testset$features <- paste(testset$epm, testset$mot_strand, testset$region, testset$type.y, sep = "_")
selected_columns <- testset[, c("features", "loc_ID", "dist_transc_border", "pred_class", "expr_class")]
selected_columns <- selected_columns[!duplicated(selected_columns),]
selected_columns02 <- testset[, c("loc_ID", "pred_class", "expr_class")]
selected_columns02 <- selected_columns02[!duplicated(selected_columns02),]
cont_table <- xtabs(~ loc_ID + features, data = selected_columns)
cont_table_df <- as.data.frame.matrix(cont_table > 0)
cont_table_df$loc_ID <- rownames(cont_table_df)
cont_table_df0 <- merge(selected_columns02, cont_table_df, by ="loc_ID")
cont_table_df_expr <- cont_table_df0[, -which(names(cont_table_df0) == "pred_class")]
colnames(cont_table_df_expr)[2] <- "target"
split<- sample.split(cont_table_df_expr, SplitRatio = 0.8)
data_train <- subset(cont_table_df_expr, split == "TRUE") 
data_test <- subset(cont_table_df_expr, split == "FALSE") 
### ### ###
#### ### ####
# Define the parameter grid for mtry
param_grid <- expand.grid(mtry = seq(1, ncol(data_train) - 1, by = 1))
# Define the control parameters for cross-validation
ctrl <- trainControl(method = "cv", number = 5)  # 5-fold cross-validation
# Perform the grid search
# Perform hyperparameter search
#model_RF <- train(
#  x = data_train[3:ncol(data_train)],
#  y = as.factor(data_train$target),
#  method = "rf",
#  trControl = ctrl,
#  tuneGrid = param_grid
#)
# Access best model and hyperparameters
#best_rf_model <- model_RF$finalModel
#best_hyperparameters <- model_RF$bestTune
#df_grid <- train(as.factor(data_train$target), 
#                 data = data_train[3:ncol(data_train)], 
#                 method = "rf",  # Random Forest method
#                 trControl = ctrl,
#                 tuneGrid = param_grid)
model_RF <- randomForest(data_train[3:ncol(data_train)],
                         as.factor(data_train$target),
                         ntree = 500,
                         trControl = ctrl,
                         tuneGrid = param_grid)
#best_rf_model <- model_RF$finalModel
#best_hyperparameters <- model_RF$bestTune
#best_rf_model
#best_hyperparameters
predictions <- predict(model_RF, data_train[3:ncol(data_train)])
baseline_accuracy <- sum(predictions == data_train$target) / nrow(data_train)
n_shuffles <- 100  # Number of shuffles
accuracy_diffs <- numeric(n_shuffles)
for (i in 1:n_shuffles) {
  shuffled_target <- sample(data_train$target)
  shuffled_predictions <- predict(model_RF, data_train[3:ncol(data_train)])
  accuracy_diffs[i] <- sum(shuffled_predictions == shuffled_target) / nrow(data_train) - baseline_accuracy
}
mean_accuracy_difference <- mean(accuracy_diffs)
predicted_classes <- predict(model_RF, data_test[3:ncol(data_test)])
predicted_classes_factor <- factor(predicted_classes, levels = levels(data_test$target))
data_test_pred <- data.frame(data_test$loc_ID, data_test$target)
data_test_pred1 <-cbind(data_test_pred, predicted_classes)
target_counts <- table(data_test_pred1$data_test.target)
Ps <- target_counts[["1"]]
Ns <- target_counts[["0"]]
predicted_counts <- table(data_test_pred1$predicted_classes)
pPs <- predicted_counts[["1"]]
pNs <- predicted_counts[["0"]]
cont_data_test_pred1f <- table(data_test_pred1$data_test.target, data_test_pred1$predicted_classes)
testC =as.data.frame(cont_data_test_pred1f)
colnames(testC) <- c("expr_class", "pred_class", "Freq")
rownames(testC) <- c("TRUE_high", "FALSE_high", "FALSE_low", "TRUE_low")
TRUE_val <- testC["TRUE_high","Freq"]+testC["TRUE_low","Freq"]
FALSE_val <- testC["FALSE_high","Freq"]+testC["FALSE_low","Freq"]
TPR <-  TRUE_val / (TRUE_val + FALSE_val)
FDR <- FALSE_val/ (FALSE_val+TRUE_val)
imp_mdl<- importance(model_RF)
imp_mdl0 <- data.frame(imp_mdl)
imp_mdl0$features<- rownames(imp_mdl0)
rownames(imp_mdl0) <- NULL
imp_mdl1 <- separate(imp_mdl0, features, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7"), sep = "_", remove = FALSE)
imp_mdl1$epm <- paste(imp_mdl1$col1, imp_mdl1$col2, imp_mdl1$col3, imp_mdl1$col4, sep = "_")
imp_mdl1$col8 <- substr(imp_mdl1$epm, nchar(imp_mdl1$epm), nchar(imp_mdl1$epm))
imp_mdl1$epm <- substr(imp_mdl1$epm, 1, nchar(imp_mdl1$epm) - 1)
imp_mdl2 <- imp_mdl1[, c("epm","col8", "col5", "col6", "col7",  "MeanDecreaseGini")]
colnames(imp_mdl2) <- c("epm", "orientation", "strand", "region", "annotation", "MeanDecreaseGini")
head(imp_mdl2)
pivot_table <- pivot_wider(imp_mdl2, id_cols = epm, names_from = c("orientation", "strand", "region", "annotation"), values_from = MeanDecreaseGini)
epm_values <- pivot_table$epm
numeric_data <- as.matrix(pivot_table[, -1])
model_RF
pdf(file = paste0(file_path_out, "RF_heatmap_output.pdf"), width = 8, height = 10)
heatmap.2(
  numeric_data,
  Rowv = TRUE,
  Colv = TRUE,
  col = (colorRampPalette(c("yellow", "purple"))(256)),
  dendrogram = "both",
  trace = "none",
  labRow = epm_values,
  cexRow = 0.7,
  cexCol = 0.6,
  main = "MeanDecreaseGini Value for EPMs as features in a RF",
  margins = c(10, 10),
  las = 2
)
dev.off()
write.table(imp_mdl0, file=paste0(file_path_out,"RF_imp_mdl.csv"), row.names = FALSE, sep = "\t",)
write.table(cont_table_df_expr, file=paste0(file_path_out,"RF_feature_space.csv"), row.names = FALSE, sep = "\t",)
output <- capture.output(print(model_RF))
file_path0 <- paste0(file_path_out,"RF_parameter_mdl.txt")
writeLines(output, con = file_path0)
tree_number <- 3  # Replace with the tree number you want (e.g., 1, 2, ...)
tree <- getTree(model_RF, k = tree_number)
length(tree)
print(head(tree))
TPR
FDR

