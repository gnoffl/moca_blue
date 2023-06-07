install.packages("AMR")
library(AMR)
##############################################################################
Filename1="TAIR10_GFF3_genes.gff _motifs_genes_GOs.txt"
##############################################################################
setwd("~/Desktop/meme_docker/meme_docker")
F1<-read.csv(Filename1,
             , sep=",", header=T, stringsAsFactors = F)
##############################################################################
colnames(F1)[13]<-"motif_type"
colnames(F1)[9]<-"strand"
F1$motif_type0<-str_extract(F1$motif_type,"_.*")
F1$motif_type1<-gsub("_","",as.character(F1$motif_type0))
F1$motif_type1<-as.numeric(F1$motif_type1)
##############################################################################
F1_TSS<-subset.data.frame(F1, F1$modistTSS== "1" & F1$plus_size=="1")
#F1_TSS<-subset.data.frame(F1, F1$modistTSS== "1")
#F1_TSS<-F1_TSS[order(F1_TSS$motif_type1),]
#F1_TSS1<-
#  F1_TSS %>%
#  filter(motif_type1%% 2==1)
#F1_TSS2<-
#  F1_TSS %>%
#  filter(motif_type1%% 2==0)
##############################################################################
F1matrix <- F1 %>% select(2,26,13,19)
F1gemo0<-table(F1matrix$GO,F1matrix$motif_type)
#head(F1gemo0)
#table(F1matrix$Modis1,F1matrix$motif_type)
ggplot(F1gemo2)
#F1gemo1<-as.data.frame(F1gemo0)
#F1gemo1<-as.table(F1gemo0)
F1gemo2<-as.data.frame.matrix(F1gemo0)
#F1gemo1<-as.matrix(F1gemo0)
g.test(F1gemo2)
chisq.test(F1gemo2$motif_1)
a<-chisq.test(F1gemo2, simulate.p.value = TRUE)

nvalues <- 3
nvars <- 2
nsamples <- 5000
data <- matrix( sample( 0:(nvalues - 1), nvars * nsamples, replace = TRUE ), nsamples, nvars )
data0<-as.data.frame.matrix(data)
g.test(data)

