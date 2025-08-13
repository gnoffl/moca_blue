library(dplyr)
# install.packages("ggplot2")
library(ggplot2)
#library(MASS)
working_directory = "~/Desktop/Rhome/moca_blue/mo_clu"
#########################################################################################################################
SEQLET_FILE = "blabla.txt"
HDF5_FILE= "Atha-0h_modisco.hdf5"
BASENAME = "SPECIES_MODEL_ETC"
#########################################################################################################################
dirpath_out = "../out/ranges"
#########################################################################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 5) {
  stop("Cannot provide more than 5 arguments!\nUsage: <script.R> [<working_directory>] [<SEQLET_FILE>] [<HDF5_FILE>] [<BASENAME>] [<dirpath_out>]", call.=FALSE)
}
if (length(args)==5) {
  dirpath_out = args[5]
}
if (length(args)>=4) {
  BASENAME = args[4]
}
if (length(args)>=3) {
  HDF5_FILE = args[3]
}
if (length(args)>=2) {
  SEQLET_FILE = args[2]
}
if (length(args)>=1) {
  working_directory = args[1]
}
# print all arguments
cat("working directory:", working_directory, "\n")
cat("SEQLET_FILE:", SEQLET_FILE, "\n")
cat("HDF5_FILE:", HDF5_FILE, "\n")
cat("BASENAME:", BASENAME, "\n")
cat("dirpath_out:", dirpath_out, "\n")

setwd(working_directory)
############################################################
#data <- read.csv("rdf5_seqlet_patternArthM0.txt"), sep=",")
data <- read.csv(file = SEQLET_FILE, sep=",")
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
  filename <- paste0("epm_", BASENAME, "_", motif, ".png")
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
result2$BASENAME <- c(BASENAME)
result2$source <- c(HDF5_FILE)

result1$BASENAME <- c(BASENAME)
result1$source <- c(HDF5_FILE)

result1$epm <- paste("epm", result1$BASENAME, result1$motif, sep="_")
result2$epm <- paste("epm", result2$BASENAME, result2$motif, sep="_")

# select only the desired columns in the result data frames
result1 <- result1 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

result2 <- result2 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

#########################################################################################################################
file_path_out <- file.path(dirpath_out, paste0(BASENAME))

write.csv(result1, file=paste0(file_path_out,"-TSS_motif_ranges_q1q9.csv"), row.names=FALSE)
write.csv(result2, file=paste0(file_path_out,"-TTS_motif_ranges_q1q9.csv"), row.names=FALSE)

# move all created images (.png files) into visualization folder
visualization_dir <- file.path(dirpath_out, "visualization")
dir.create(visualization_dir, recursive = TRUE, showWarnings = FALSE)
png_files <- list.files(pattern = "\\.png$", full.names = TRUE)
for (file in png_files) {
  file.rename(file, file.path(visualization_dir, basename(file)))
}
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
