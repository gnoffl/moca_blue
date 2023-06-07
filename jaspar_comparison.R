install.packages("BiocManager")
BiocManager::install("JASPAR2022")



library(JASPAR2022)
library(TFBSTools)
library(rhdf5)
library(seqLogo)

# Getting motifs from JASPAR
jaspar_motifs_plants <- getMatrixSet(JASPAR2022,
                                     opts = list(tax_group='plants',
                                                 matrixtype='PWM'))

jaspar_pwm <- list()
jaspar_pwm_ids <- c()
jaspar_fam_name <- c()
for (name_idx in 1:length(names(jaspar_motifs_plants))) {
  jaspar_pwm[[name_idx]] <- as.matrix(jaspar_motifs_plants[[
    names(jaspar_motifs_plants)[name_idx]]])
  jaspar_pwm_ids[name_idx] <- names(jaspar_motifs_plants)[name_idx]
  jaspar_fam_name[name_idx] <- name(jaspar_motifs_plants)[name_idx]
  
}

meta_info <- do.call(cbind, list(family=jaspar_fam_name, id=jaspar_pwm_ids))
write.csv(meta_info, 'jaspar_meta_info.csv', row.names = F)
for (file_name in list.files('modisco')) {
  pfms <- h5read(paste0('modisco/', file_name), 'pfms')
  mot_sim_list <- list()
  
  for (idx in seq(dim(pfms)[3])) {
   mot_to_jas <- c()
   pfm <- as.matrix(pfms[, , idx]) 
   pfm <- log2(((pfm/colSums(pfm)) + 0.00001)/0.25)
   rownames(pfm) <- c('A', 'C', 'G', 'T')
   for (idx_2 in seq(length(jaspar_motifs_plants))) {
     mot_to_jas[idx_2] <- PWMSimilarity(pfm, jaspar_pwm[[idx_2]], method = 'Pearson')
     
   }
   mot_sim_list[[idx]] <- mot_to_jas
  }
  sim_mat<-do.call(rbind, mot_sim_list)
  colnames(sim_mat) <- jaspar_pwm_ids
  write.table(sim_mat, gsub('pred_motifs.h5', 'sim.csv', file_name),
            row.names = F, sep = '\t')
  
}


