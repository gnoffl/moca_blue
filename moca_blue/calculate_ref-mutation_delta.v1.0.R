setwd("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/")
# # # # # # # # # # # # # # # # #
library(dbplyr)
library(ggplot2)

 # # # # # # # # # # # # # # # # 
df_ref_list <- c("Atha-24h_logmaxTPM.2.csv")
#
df_predict_list <- c("Atha-24h_model__Atha_predict.csv",
             "Atha-24h_model__AlyrSNPs_predict.csv",
             "Atha-24h_model__0.03SNPs_predict.csv",
             "Atha-24h_model__0.045SNPs_predict.csv",
             "Atha-24h_model__0.06SNPs_predict.csv")
df_names <- c("Atha_predict",
              "AlyrSNPs_predict",
              "0.03SNPs_predict",
              "0.045SNPs_predict",
              "0.06SNPs_predict")
# # # # # # IMPORT # # # # # # # # 
for (i in 1:length(df_ref_list)) {
  assign(gsub(".csv", "", df_ref_list[i]), read.csv(df_ref_list[i]))
}
# # # # # # # # # # # # # # # # # # 
for (i in 1:length(df_predict_list)) {
  assign(gsub(".csv", "", df_predict_list[i]), read.csv(df_predict_list[i]))
}
# # # # # # WRANGLE # # # # # # # #
ref_df01<- `Atha-24h_logmaxTPM.2`
colnames(ref_df01)[1]<- "gene_ids"
lower_quantile <- quantile(ref_df01$logMaxTPM, 0.25)
upper_quantile <- quantile(ref_df01$logMaxTPM, 0.75)
ref_df01$expr_class <- ifelse(ref_df01$logMaxTPM > upper_quantile, "1", 
                              ifelse(ref_df01$logMaxTPM > lower_quantile, "2","0"))
ref_df01<- ref_df01 %>%
  select(gene_ids, logMaxTPM, expr_class)
# # # # # # # # # # # # # # # # # # #
for (i in 1:length(df_predict_list)) {
  df <- read.csv(df_predict_list[i])
  colnames(df) <- c("gene_ids", df_names[i]) # Rename the second column to match the corresponding df name
  assign(gsub(".csv", "", df_predict_list[i]), df) # Assign the renamed df back to the environment
}

R_mergdf01 <- ref_df01 # Start with the first data frame
for (df_name in df_predict_list) {
  df <- get(gsub(".csv", "", df_name)) # Get the data frame by name
  R_mergdf01 <- merge(R_mergdf01, df, by = "gene_ids", all = TRUE) # Merge with previous merged_df
}
# # # # # # # # # # # # # # # # # # # #
# Calculate variances
R_mergdf01$delta_Atha_Atha <- (R_mergdf01$Atha_predict - R_mergdf01$Atha_predict)*-1
R_mergdf01$delta_Atha_Alyr <- (R_mergdf01$Atha_predict - R_mergdf01$AlyrSNPs_predict)*-1
R_mergdf01$delta_Atha_0.03SNPs <- (R_mergdf01$Atha_predict - R_mergdf01$`0.03SNPs_predict`)*-1
R_mergdf01$delta_Atha_0.045SNPs <- (R_mergdf01$Atha_predict - R_mergdf01$`0.045SNPs_predict`)*-1
R_mergdf01$delta_Atha_0.06SNPs <- (R_mergdf01$Atha_predict - R_mergdf01$`0.06SNPs_predict`)*-1
#
averages <- aggregate(. ~ R_mergdf01$expr_class, data = R_mergdf01[, grepl("^delta_", names(R_mergdf01))], mean)
# # # # # # # # ##
# Calculate summary statistics for each column starting with "delta_"
stats <- sapply(R_mergdf01[, grepl("^delta_", names(R_mergdf01))], summary)
#stats_mean_sd <- sapply(R_mergdf01[, grepl("^delta_", names(R_mergdf01))], function(x) c(mean = mean(x), sd = sd(x)))
# Convert the result to a data frame
stats_df <- as.data.frame(stats)
# Transpose the data frame for better readability
stats_df <- as.data.frame(t(stats_df))
#################
# Reshape data for plotting
delta_data_long <- reshape2::melt(R_mergdf01, id.vars = "expr_class", measure.vars = grep("^delta_", names(R_mergdf01), value = TRUE))
# Create separate violin plots for each column starting with "delta_"
p <- ggplot(delta_data_long, aes(x = expr_class, y = value, fill = expr_class)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "expr_class", y = "Delta(*-1)", title = "Distribution of delta values by expr_class") +
  facet_grid(. ~ variable, scales = "free_y", switch = "x") +  # Arrange in a row
  theme_minimal() +
  theme(strip.text.x = element_text(angle = 45, hjust = 1))

# Add summary statistics
p0<- p + stat_summary(fun.data = function(x) {
  y <- mean(x)
  ymin <- y - sd(x)
  ymax <- y + sd(x)
  return(data.frame(y = y, ymin = ymin, ymax = ymax))
}, geom = "crossbar", width = 0.5, color = "black", alpha = 0.5)
##################
ggsave("Atha-24h_ref-mutation_delta_distrib.png", p0, width = 12, height = 4, units = "in")
write.csv(averages, "Atha-24h_ref-mutation_delta_averages_predict.csv" )
write.csv(stats_df, "Atha-24h_ref-mutation_delta_sumstats.csv")

head(AthAly.genome.snp)
data_df <- as.data.frame(AthAly.genome.snp)
gene_counts <- as.data.frame(table(data_df$gene))
sum(gene_counts$Freq)/nrow(gene_counts)
count_equal <- sum(data_df$tha == data_df$lyr)
mutations<-nrow(data_df)-count_equal
#number of genes = 27655
#average gene length = 2373.544
av_mut_per_gene<-mutations/27655
av_mut_per_gene/2373.544
#avergae_mutation_per_gene 0.01939919