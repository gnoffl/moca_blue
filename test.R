library(TFBSTools)

FILE1 = "all_motifs_20230505.jaspar"
#######################################################
dirpath_in = "../Mo_Nom/out/"
dirpath_out = "./out/"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
File1 <- paste0(dirpath_in,FILE1)
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################################################
cwm1 <- read_jaspar(File1)
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  
  # Extract the motif name and the number after the last "_" underscore
  motif_name <- attr(cwm1[[i]], "name")
  nsites <- as.numeric(sub(".+_(\\d+)$", "\\1", motif_name))
  
  # Assign the nsites value to the "nsites" field of the motif object
  cwm1[[i]]["nsites"] <- nsites
}
##################################################
# Loop through each motif object in the list
for (i in seq_along(cwm1)) {
  
  # Extract the Total IC and the Consensus values from the motif object
  total_ic <- attr(cwm1[[i]], "icscore")
  total_ic_rounded <- round(total_ic, 1)
  
  consensus <- attr(cwm1[[i]], "consensus")
  
  # Combine the Total IC and the Consensus values separated by "_" to the motif name
  motif_name <- attr(cwm1[[i]], "name")
  new_motif_name <- paste0(motif_name, "_", total_ic_rounded, "_", consensus)
  
  # Assign the new motif name to the motif object
  attr(cwm1[[i]], "name") <- new_motif_name
}
##################################################
pwm_uni0<-convert_motifs(
  cwm1, class = "TFBSTools-PWMatrix")

pcm<-convert_motifs(
  cwm1, class = "motifStack-pcm")
##################################################
sum<-as.data.frame(summarise_motifs(pcm))
write.csv(sum, file = paste0(dirpath_out,FILE1,"summary-test.txt"))
##################################################
#c_pcm<-clusterMotifs(pcm) ### !!! TIME TO GET A COFFEEE !!! ###
#hc<- c_pcm
#motifs<-pcm[hc$order]
##################################################
#write.tree(as.phylo(c_pcm), file = paste0(dirpath_out,FILE1,"-SW.nwk"))

########################################################
#motifs<-pcm[hc$order]
#motifs <- lapply(motifs, pcm2pfm)
#d1o alignment
#compare_motif()
##########################################################
cwm1[[1]]
##########################################################
c<-compare_motifs(cwm1, method = "EUCL")
c0<-as.data.frame(c)
write.csv(c0, file = paste0(dirpath_out,FILE1,"_matrix-EUCL-test.txt"))
#assign scores from data.frame to branches for selection
tree <- motif_tree(cwm1, layout = "rectangular", db.scores = "scores", method = "EUCL")
###############################################################
hc <- clusterMotifs(pcm)
phylog <- hclust2phylog(hc)
leaves <- names(phylog$leaves)
pfms <- pfms[leaves]
## create a list of pfm objects
pfms <- mapply(pfms, names(pfms),
               FUN=function(.pfm, .name){
                 new("pfm",mat=.pfm, name=.name)})
## extract the motif signatures
motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.0001, min.freq=1)
sig <- signatures(motifSig)
## get the group color for each signature
gpCol <- sigColor(motifSig)

#library(RColorBrewer)
#color <- brewer.pal(12, "Set3")
## plot the logo stack with pile style.
#motifPiles(phylog=phylog, pfms=pfms, pfms2=sig,
#           col.tree=rep(color, each=5),
#           col.leaves=rep(rev(color), each=5),
#           col.pfms2=gpCol,
#           r.anno=c(0.02, 0.03, 0.04),
#           col.anno=list(sample(colors(), 50),
#                         sample(colors(), 50),
#                         sample(colors(), 50)),
#           motifScale="logarithmic",
#           plotIndex=TRUE)