library(dplyr)
#library(MASS)

#########################################################################################################################
NAME0="Rdf5_seqlets_pattern"
SPEC="Soly"
MODEL="MSR"
FILE= "Solanum_MSR_modisco.hdf5"
#########################################################################################################################
data <- read.csv(file = paste0(NAME0,SPEC,MODEL,".txt"), sep=",")
#########################################################################################################################
########################################################################################################################
#########################################################################################################################
data$motif <- ifelse(data$metacluster == "metacluster_0", "p0", "p1")
data$motif <- paste0(data$motif, "m", sprintf("%02d", as.numeric(substring(data$pattern, 9))))
#########################################################################################################################

#data$start <- as.numeric(data$start)
range1 <- data %>% filter(start >= 1 & end <= 1500)
range2 <- data %>% filter(start >= 1520 & end <= 3000)

range1$trunc_start <- floor(range1$start / 10)
range2$trunc_start <- floor(range2$start / 10)
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
            q10 = quantile(start, 0.10),
            median = median(start),
            q90 = quantile(start, 0.90),
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
            q10 = quantile(start, 0.10),
            median = median(start),
            q90 = quantile(start, 0.90),
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

result1$epm <- paste("epm", result1$Species, substr(result1$Model, 1, 1), result1$motif, sep="_")
result2$epm <- paste("epm", result2$Species, substr(result2$Model, 1, 1), result2$motif, sep="_")

# select only the desired columns in the result data frames
result1 <- result1 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

result2 <- result2 %>%
  select(epm, min, max, mean, median, mode, q10,  q90,  sd, cv, iqr, number, source)

#########################################################################################################################

write.csv(result1, file=paste0(SPEC,MODEL,"-TSS_motif_position.csv"), row.names=FALSE)
write.csv(result2, file=paste0(SPEC,MODEL,"-TTS_motif_position.csv"), row.names=FALSE)