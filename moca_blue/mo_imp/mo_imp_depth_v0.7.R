######################
library(rhdf5)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
library(dbscan)
library(plotly)
library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(circlize)
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
model_parameter$abs_sum_imp_score <- rowSums(abs(imp_scores))
model_parameter$TRUE_pred_imp <- ifelse(model_parameter$sum_imp_score2 == model_parameter$true_target, TRUE, FALSE)
gene_imp_scores <- as.data.frame(model_parameter$gene_id)
gene_imp_scores <- cbind(gene_imp_scores, imp_scores)
colnames(gene_imp_scores)[1] <- "loc_ID"
########################################### Calculate and visualize the correlation between predicted probability and imp_scores
model_parameter_0 <- model_parameter[, c(1,2,4,7,9)] # select values of interest
plot(model_parameter_0$pred_prob, model_parameter_0$sum_imp_score,
     pch = 16, cex = 0.8, xlab = "pred_prob", ylab = "sum_imp_score", main = "Scatter Plot", col = "black")
# Assuming model_parameter_0 is your data frame
plot(model_parameter_0$pred_prob, model_parameter_0$abs_sum_imp_score,
     pch = 16, cex = 0.8, xlab = "pred_prob", ylab = "abs_sum_imp_score", main = "Scatter Plot", col = "black")
# Add line for the median of abs_sum_imp_score
abline(h = median(model_parameter_0$abs_sum_imp_score), col = "red", lty = 3)
# Add horizontal line at y = 0.5
abline(h = 1, col = "blue", lty = 2)
#############################################
data <- cbind(model_parameter_0$pred_prob, model_parameter_0$abs_sum_imp_score)
# Standardize the data (important for DBSCAN)
scaled_data <- scale(data)
# Combine the standardized features into a single matrix
features_matrix <- as.matrix(scaled_data)
# Applying DBSCAN
epsilon <- 0.1  # The maximum distance between two samples for one to be considered as in the neighborhood of the other
min_samples <- 100  # The number of samples (or total weight) in a neighborhood for a point to be considered as a core point
dbscan_result <- dbscan(features_matrix, eps = epsilon, MinPts = min_samples)
# Visualizing the results
plot(features_matrix, col = dbscan_result$cluster + 1, pch = 16, main = "DBSCAN Clustering", xlab = "pred_prob", ylab = "abs_sum_imp_score")
legend("topright", legend = unique(dbscan_result$cluster), col = unique(dbscan_result$cluster) + 1, pch = 16, title = "Cluster")
# Filter data based on cluster assignments
cluster_to_remove <- 0
filtered_data <- as.data.frame(data[dbscan_result$cluster != cluster_to_remove, ])
filtered_data$group <- ifelse(filtered_data$V1 > 0.5, "1", "0")
head(filtered_data)
#filtered_data$group <- factor(filtered_data$group, levels = c("Group1", "Group2"))
boxplot(V2 ~ group, data = filtered_data, col = c("red", "lightgreen"), main = "Boxplot for Groups 0 and 1", ylab = "V2")
#legend("topright", legend = levels(filtered_data$group), fill = c("red", "lightgreen"), title = "Groups")

dbscan_table <- data.frame(Index = seq_len(nrow(features_matrix)),
                           DBSCluster = dbscan_result$cluster)

model_parameter_0$Index <- seq_len(nrow(model_parameter_0)) # select values of interest
model_parameter_1 <- merge(model_parameter_0,dbscan_table, by= "Index")
############################################
DfC_high <- model_parameter_1[model_parameter_1$DBSCluster == 2, ]
Dfc_high_scores <- imp_scores[DfC_high$Index, ]
Dfc_high_scores0 <- abs(Dfc_high_scores)
#############
Dfc_high_scores <- as.matrix(Dfc_high_scores)

cor_matrix <- cor(t(Dfc_high_scores))

order_indices <- order(rowSums(1 - cor_matrix))
Dfc_high_scores_order <- Dfc_high_scores[order_indices, ]


# Perform hierarchical clustering
clustering_result <- hclust(as.dist(1 - cor_matrix), method = "complete")
###
par(mfrow=c(1,1))
# Plot the dendrogram
plot(clustering_result, main = "Hierarchical Clustering", xlab = "Arrays", sub = "Dendrogram")
#
#heatmap(
#  Dfc_high_scores_order,
#  col = circlize::colorRamp2(c(min(Dfc_high_scores_order), 0, max(Dfc_high_scores_order)), c("cyan", "white", "darkorange")),
#  cluster_rows = FALSE, cluster_columns = FALSE,
#  show_heatmap_legend = FALSE,
#  use_raster = TRUE,
#  raster_resize = TRUE
#)

heatmap(
  Dfc_high_scores,
#  col = circlize::colorRamp2(c(min(Dfc_high_scores_order), 0, max(Dfc_high_scores_order)), c("cyan", "white", "darkorange")),
  Rowv = as.dendrogram(hclust(dist(Dfc_high_scores))),
  Colv = NA,
  use_raster = TRUE,
  raster_resize = TRUE,
main = "Heatmap with Hierarchical Clustering", xlab = "Arrays", ylab = "Genes")
########################
# Volcano dataset
#volcano

# Heatmap 
Dfc_high_scores_order %>%
  as_tibble() %>%
  rowid_to_column(var = "X") %>%
  gather(key = "Y", value = "Z", -1) %>%
  mutate(Y = as.numeric(gsub("V", "", Y))) %>%
  
  # Order rows based on similarity
  arrange(X) %>%
  
  ggplot(aes(X, Y, fill = Z)) +
  geom_tile() +
  theme_minimal() +  # Use your desired theme
  theme(legend.position = "none")
#########################

num_clusters <- 3 # You can choose the number of clusters
clusters <- cutree(clustering_result, num_clusters)

# Group data based on clusters
grouped_data <- split(Dfc_high_scores, clusters)

# Average values within each cluster
averaged_data <- lapply(grouped_data, function(cluster_data) colMeans(cluster_data))

# Convert the list of averages back to a matrix
averaged_matrix <- do.call(cbind, averaged_data)

# Normalize the averaged values (similar to your example)
column_sum <- apply(averaged_matrix, 2, sum)
normalized_averages <- sapply(averaged_matrix, function(col) col / abs(column_sum))

# Plot the normalized averaged values
# Calculate the y-axis ticks at 500 intervals
x_ticks <- seq(1, ncol(normalized_averages), by = 500)

# Plot the normalized averaged values with custom y-axis ticks
par(mfrow = c(2, 1))  # Adjust the layout for a clearer plot
matplot(t(normalized_averages), type = "l", col = 1:num_clusters,
        xlab = "Column Index", ylab = "Normalized Average Value",
        main = "Normalized Averaged Values within Clusters", ylim = c(0, max(normalized_averages)),
        yaxt = "n")  # Suppress the y-axis

# Add custom y-axis ticks
axis(1, at = x_ticks, labels = x_ticks)
# Add legend
#legend("topright", legend = paste("Cluster", 1:num_clusters), col = 1:num_clusters, lty = 1)
##########################

head(DfC_high)

############################################

selected_gene <- imp_scores[1474,]
slctd_gene_t0 <- as.data.frame(t(selected_gene)) #Change here
column_1_sum <- sum(slctd_gene_t0[, 1])
column_1_sum
rooted_sqrd_val <- sqrt(slctd_gene_t0[, 1]^2)
rtd_sqr_sum <- sum(rooted_sqrd_val)
rtd_sqr_sum
abs_imp_sum <- sum(abs(slctd_gene_t0[, 1]))
#abs_imp_sum
slctd_gene_t0$norm_imp <- (slctd_gene_t0[, 1] / abs(column_1_sum))
slctd_gene_t0$norm_prd_imp <- c(sqrt(slctd_gene_t0[, 1]^2)/abs_imp_sum)
head(slctd_gene_t0)
plot(1:nrow(slctd_gene_t0), slctd_gene_t0[, 1], type = "l", col = rgb(0, 0, 1, alpha = 0.3), xlab = "Row", ylab = "Values", main = "Line Plot of Columns")
lines(1:nrow(slctd_gene_t0), slctd_gene_t0$norm_imp, col = rgb(0, 1, 0, alpha = 0.3))
lines(1:nrow(slctd_gene_t0), slctd_gene_t0$norm_prd_imp, col = rgb(0, 0, 0, alpha = 0.8))

#legend("topright", legend = c("1", "norm_imp", "norm_prd_imp"), col = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5), rgb(0, 1, 0, alpha = 0.5)), lty = 1)
######################################



#######################################
# Selecting the third row and transposing it to ensure a numeric vector
#selected_gene_fft <- fft(as.numeric(t(imp_scores[1444, ])))

# Extracting magnitude information
#magnitude <- Mod(selected_gene_fft)

# Plotting the magnitude spectrum
#plot(magnitude, type = "l", xlab = "Frequency", ylab = "Magnitude", main = "Magnitude Spectrum")

