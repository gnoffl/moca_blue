# Set the working directory
#setwd("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/")
# List all files with "results" in their name and sort them alphanumerically
file_list <- sort(list.files(pattern = "result2"))

# Create a list to store data frames
data_list <- list()

# Import data from each file and store in the list
for (file in file_list) {
  data <- read.csv(file)
    data <- data[data$val_acc >= 0.5, ]
  
  data_list[[file]] <- data
}
num_files <- length(file_list)
num_rows <- ifelse(num_files > 4, ceiling(num_files / 2), num_files)
num_cols <- ifelse(num_files > 4, 2, 1)
par(mfrow=c(num_rows, num_cols), mar=c(2, 2, 2, 2), oma=c(0, 0, 2, 0))
for (i in 1:length(file_list)) {
  file <- file_list[i]
  data <- data_list[[file]]
  colors <- colorRampPalette(c("lightblue", "lightgrey"))(nrow(data$val_acc))
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0.5, length(data$val_acc) + 0.5), ylim = c(0.5, 1), main = file)
  rect(xleft = seq(0.5, length(data$val_acc) - 0.5), ybottom = 0.5, xright = seq(1.5, length(data$val_acc) + 0.5), ytop = data$val_acc, col = colors, border = NA)
  axis(1, at = seq(1, length(data$val_acc), 1), labels = FALSE)
  axis(2)
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2), oma=c(0, 0, 0, 0))
