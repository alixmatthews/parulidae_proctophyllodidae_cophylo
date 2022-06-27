#### R CODE FOR PAIRWISE GENETIC DISTANCES ####


#### SET WORKING DIRECTORY
setwd("/vol_c/")


#### LOAD LIBRARIES 
library(ape)
library(phangorn)

#### distances
# COI
COI <- read.dna(file="/vol_c/mt_genes/20211115/aligned_fasta/concat/coi_cat.fasta.mafftaligned.fasta", "f")

# -- P-dist
COI.Pdist<-dist.dna(COI, "raw", pairwise.deletion = TRUE, as.matrix=TRUE) 
View(COI.Pdist)
write.table(COI.Pdist, file = "COI.Pdist.csv", sep=",")

# -- K2P dist
COI.K2Pdist<-dist.dna(COI, "K80", pairwise.deletion = TRUE, as.matrix=TRUE) 
View(COI.K2Pdist)
write.table(COI.K2Pdist, file = "COI.K2Pdist.csv", sep=",")

# -- JC dist
COI.JCdist<-dist.dna(COI, "JC69", pairwise.deletion = TRUE, as.matrix=TRUE) 
write.table(COI.JCdist, file = "COI.JCdist.csv", sep=",")




