setwd("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/")
# # # # # # # # # # # # # # # # #
library(dbplyr)
library(ggplot2)

 # # # # # # # # # # # # # # # # 
df_ref_list <- c("Atha-12h_logmaxTPM.2.csv")
#
df_predict_list <- c("Atha-12h_model__Atha_predict.csv",
             "Atha-12h_model__AlyrSNPs_predict.csv",
             "Atha-12h_model__0.03SNPs_predict.csv",
             "Atha-12h_model__0.045SNPs_predict.csv",
             "Atha-12h_model__0.06SNPs_predict.csv")
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
ref_df01<- `Atha-12h_logmaxTPM.2`
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
#######################################

R_mergdf02 <- R_mergdf01[R_mergdf01$expr_class != 2, ]
R_mergdf02$TF_Atha_predict <- ifelse((R_mergdf02$expr_class == 1 & R_mergdf02$Atha_predict > 0.5) |
                                              (R_mergdf02$expr_class == 0 & R_mergdf02$Atha_predict < 0.5), "TRUE", "FALSE")
R_mergdf02$TF_AlyrSNPs_predict <- ifelse((R_mergdf02$expr_class == 1 & R_mergdf02$AlyrSNPs_predict > 0.5) |
                                            (R_mergdf02$expr_class == 0 & R_mergdf02$AlyrSNPs_predict < 0.5), "TRUE", "FALSE")
R_mergdf02$TF_0.03SNPs_predict <- ifelse((R_mergdf02$expr_class == 1 & R_mergdf02$'0.03SNPs_predict' > 0.5) |
                                       (R_mergdf02$expr_class == 0 & R_mergdf02$'0.03SNPs_predict' < 0.5), "TRUE", "FALSE")
R_mergdf02$TF_0.045SNPs_predict <- ifelse((R_mergdf02$expr_class == 1 & R_mergdf02$'0.045SNPs_predict' > 0.5) |
                                           (R_mergdf02$expr_class == 0 & R_mergdf02$'0.045SNPs_predict' < 0.5), "TRUE", "FALSE")
R_mergdf02$TF_0.06SNPs_predict <- ifelse((R_mergdf02$expr_class == 1 & R_mergdf02$'0.06SNPs_predict' > 0.5) |
                                            (R_mergdf02$expr_class == 0 & R_mergdf02$'0.06SNPs_predict' < 0.5), "TRUE", "FALSE")
#######################################
# Select columns starting with "TF_"
tf_columns <- grep("^TF_", names(R_mergdf02), value = TRUE)

# Function to calculate percentages and FDR for a single column
calculate_statistics <- function(column_name) {
  # Check if column contains only NA values
  if(all(is.na(R_mergdf02[[column_name]]))) {
    return(data.frame(
      Column_Name = column_name,
      Percentage_TRUE = NA,
      Percentage_FALSE = NA,
      FDR = NA
    ))
  }
  
  # Calculate percentage of TRUE and FALSE values
  percentage_true <- mean(R_mergdf02[[column_name]] == "TRUE", na.rm = TRUE) * 100
  percentage_false <- 100 - percentage_true
  # Calculate False Discovery Rate (FDR)
  fdr <- percentage_false / (percentage_true + percentage_false)
  return(data.frame(
    Column_Name = column_name,
    Percentage_TRUE = percentage_true,
    Percentage_FALSE = percentage_false,
    FDR = fdr
  ))
}
# Calculate statistics for each column separately
results <- lapply(tf_columns, calculate_statistics)
# Combine results into a single dataframe
results_df <- do.call(rbind, results)
# Print the results
# Export results_df as a CSV file
write.csv(results_df, "Atha-12h_ref-mut-TDRs_stats.csv", row.names = FALSE)
# Export R_mergdf02 as a CSV file
write.csv(R_mergdf02, "Atha-12h_ref-mut-TDRs.csv", row.names = FALSE)

