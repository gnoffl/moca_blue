library(dplyr)
library(ggplot2)
#library(MASS)
setwd("~/Desktop/Rhome/moca_blue/mo_clu")
#########################################################################################################################
NAME0="rdf5_seqlet_pattern"
SPEC="Arth"
MODEL="D0"
FILE= "Atha-0h_modisco.hdf5"
#########################################################################################################################
dirpath_in = "./out"
dirpath_out = "./out"
############################################################
file_path_in <- file.path(dirpath_in, paste0(NAME0,SPEC,MODEL))
#data <- read.csv("rdf5_seqlet_patternArthM0.txt"), sep=",")
data <- read.csv(file = paste0(file_path_in,".txt"), sep=",")
#########################################################################################################################
########################################################################################################################
#########################################################################################################################
data$motif <- ifelse(data$metacluster == "metacluster_0", "p0", "p1")
data$motif <- paste0(data$motif, "m", sprintf("%02d", as.numeric(substring(data$pattern, 9))))
#################################################################################
####
unique_motifs <- unique(data$motif)

# Iterate through each unique motif
for (motif in unique_motifs) {
  motif_data <- data[data$motif == motif, ]  # Subset data for current motif
  total_motifs <- nrow(motif_data)
  
  # Create ggplot
  p <- ggplot(motif_data, aes(x = start)) +
    geom_histogram(binwidth = 20, aes(y = ..count.. / total_motifs * 100), fill = "darkgrey", color = "white") +
    labs(title = paste("EPM Occurrence for Motif:", motif),
         x = "Position",
         y = "Percentage of Occurrence") +
    scale_y_continuous(breaks = c(0, 5, 10, 20, 25))  +
    theme(panel.grid.major = element_line(color = "lightgrey"),
          panel.background = element_rect(fill = "white"))
  
  # Set x-axis limits
  p <- p + xlim(0, 3020)
  #  + scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2020, 2520, 3020))
  
  # Save the plot
  filename <- paste0("epm_", SPEC, "_", MODEL, "_", motif, ".png")
  ggsave(filename, p, width = 10, height = 2.5, units = "in", dpi = 300)
}

#########################################################################################################################
#data$start <- as.numeric(data$start)
range1 <- data %>% filter(start >= 1 & end <= 1500)
range2 <- data %>% filter(start >= 1520 & end <= 3000)
range1$trunc_start <- floor(range1$start / 10)
range2$trunc_start <- floor(range2$start / 10)
range2$start <- c(3020-range2$start)
range2$end <- c(3020-range2$end)
#####################################################
##################################################### Create function to calculate mode #################################

customMode <- function(x) {
  freq <- table(x)
  mode <- as.numeric(names(freq)[which.max(freq)])
  return(mode)
}
#########################################################################################################################

# calculate the summary statistics, including the mode and its frequency
result1 <- range1 %>%
  group_by(motif) %>%
  summarize(min = min(start),
            max = max(start),
            q10 = quantile(start, 0.1),
            median = median(start),
            q90 = quantile(start, 0.9),
            mode = customMode(trunc_start),
            mean = mean(start),
            sd = sd(start),
            cv = sd(start) / mean(start) * 100,
            iqr = q90 - q10,
            number = n())

# extract the mode value and frequency from the mode vector and multiply the mode by 10
result1$mode <- c(result1$mode * 10)  #mode uses decimal

result2 <- range2 %>%
  group_by(motif) %>%
  summarize(min = min(start),
            max = max(start),
            q10 = quantile(start, 0.1),
            median = median(start),
            q90 = quantile(start, 0.9),
            mode = customMode(trunc_start),
            mean = mean(start),
            sd = sd(start),
            cv = sd(start) / mean(start) * 100,
            iqr = q90 - q10,
            number = n())

# extract the mode value and frequency from the mode vector and multiply the mode by 10
result2$mode <- c(result2$mode * 10)  #mode uses decimal

# add additional information to the result data frames
result2$Species <- c(SPEC)
result2$Model <- c(MODEL)
result2$source <- c(FILE)

result1$Species <- c(SPEC)
result1$Model <- c(MODEL)
result1$source <- c(FILE)

result1$epm <- paste("epm", result1$Species, result1$Model, result1$motif, sep="_")
result2$epm <- paste("epm", result2$Species, result2$Model, result2$motif, sep="_")

# select only the desired columns in the result data frames
result1 <- result1 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

result2 <- result2 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

#########################################################################################################################
file_path_out <- file.path(dirpath_out, paste0(SPEC,MODEL))

write.csv(result1, file=paste0(file_path_out,"-TSS_motif_ranges_q1q9.csv"), row.names=FALSE)
write.csv(result2, file=paste0(file_path_out,"-TTS_motif_ranges_q1q9.csv"), row.names=FALSE)
#####################
#motif_data <- subset(data, motif == "p0m00")
#total_motifs <- nrow(motif_data)
#p <- ggplot(motif_data, aes(x = start)) +
#  geom_histogram(binwidth = 20, aes(y = ..count.. / total_motifs * 100), fill = rgb(0.218,0.7,0.59), color = "white") +
#  labs(title = "epm occurence",
#       x = "position",
#       y = "percentage of occ.") +
#  xlim(0, 3020)
#filename <- paste0("epm_", SPEC, "_", MODEL, "_", motif_data$motif[1], ".png")
#ggsave(filename, p, width = 8, height = 6, units = "in", dpi = 300)
