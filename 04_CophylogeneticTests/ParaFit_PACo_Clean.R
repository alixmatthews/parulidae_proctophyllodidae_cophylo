#### ParaFit and PACo analyses
#### Alix Matthews
#### Initially done in Fall 2021, cleaned 28 June 2022
#### I cleaned the code so that results are reproducible and easy to follow
#### R version 3.6.3


#### INITIALIZATION ####

# Set wd
setwd("/Users/alix22792/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_phylo/Jetstream")

# Load libraries
library(ape) # version 5.5
library(phangorn) # version 2.7.0
library(parallel) # version 3.6.3
library(vegan) # version 2.5-7
library(paco) # version 0.4.2
# Phytools downloaded from GitHub (not CRAN): https://github.com/liamrevell/phytools
library(phytools) # version 0.7-90

# set alpha level
p_level<-0.0549







#### LOAD AND CLEAN UP TREES ####

#### ~~ MITE TREE - CONCATENATED (AICc) ####
# Read in 'parasite' tree (rooted)
original.p.tree<-read.nexus("/Users/alix22792/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_phylo/Jetstream/mt_genes/20211115/aligned_fasta/concat/Concatenated_AICc/Mites_tree_concat_AICc_rooted.treefile") # this is a nexus format. Otherwise, use "read.tree"

# view it to check it
plot(original.p.tree)

# remove 'parasite' outgroups
no.outgroup.p.tree <- drop.tip(original.p.tree, c('YEWA_941_P_Pr_CTATCCAC-TGTGCGTT', 'YEWA_941_P_Tr_ACGCTTCT-GACTTAGG', 'Proctophyllodes_miliariae'))

# view it to check it
plot(no.outgroup.p.tree)

# Trim tree to one tip per OTU as specified by species delineation (in this case, ABGD)
# first, read in tip labels
no.outgroup.p.tree$tip.label

# drop tips (randomly per OTU)
ABGD.OTU.p.tree <- drop.tip(no.outgroup.p.tree, c('BTNW_765_P_CCAAGTAG-TCACGTTC','YTWA_970_R_GAACCTTC-GTACCTTG','CMWA_952_R_CTAAGACC-CTGTTAGG','CSWA_939_R_TCCATTGC-TCCTACCT','BWWA_412_R_GTACACCT-TTGGTCTC',"BLBW_861_R_AACGCCTT-TTCTCTCG","GWWA_928_R_ACGAGAAC-CGAACTGT","CAWA_366_A_CCGGAATA-CCTACTGA_dd","BTNW_492_R_AAGGACCA-GTGTTCCT","BTBW_850_R_AGGTTCCT-AGGTTCGA","AMRE_967_R_ACGGACTT-TCTCCGAT","CERW_358_A_TCTTCGAC-GAGCAGTA","BPWA_337_C_CAAGCCAA-CTTCGTTC_dd","WEWA_727_R_ACCATCCT-GTTCATGG","PIWA_948_R_ATAACGCC-TTGGTGAG","MOWA_611_R_HC7TNCCX2_subsampled_dd","COYE_868_R_CCAACACT-AGTTCGTC","NAWA_328_A_CCACAACA-CTTGGATG","KEWA_584_R_TCTAGGAG-CGTTGCAA","LOWA_438_A_AACTTGCC-GAAGTTGG"))

# view it
plot(ABGD.OTU.p.tree)


# Change tip labels on this OTU tree
# first, read the remaining tips in to see the correct order needed for relabeling
ABGD.OTU.p.tree$tip.label

# relabel them
new_tiplabels.p.tree <- c("Tyrannidectes_charitomenos", "Tyrannidectes_aff_charitomenos", "Amerodectes_ischyros", "Amerodectes_aff_ischyros", "Amerodectes_helmitheros", "Amerodectes_sp_n_1", "Amerodectes_hribari", "Amerodectes_protonotaria", "Amerodectes_seiurus", "Amerodectes_jonesborensis")

# assign new labels to tree
ABGD.OTU.p.tree$tip.label <- new_tiplabels.p.tree

# plot it to check it
plot(ABGD.OTU.p.tree)

#### Let's export this tree for future use (eg in JANE, empress)
ape::write.tree(ABGD.OTU.p.tree, file = "ABGD.OTU.p.tree.txt")
ape::write.nexus(ABGD.OTU.p.tree, file = "Mites_concat_AICc.nex")
ape::write.tree(ABGD.OTU.p.tree, file = "Mites_concat_AICc.nwk")



#### ~~ MITE TREE - COI ONLY ####
#### Read in 'parasite' tree (rooted)
original.coi.p.tree<-read.nexus("/Users/alix22792/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_phylo/Jetstream/mt_genes/20211115/aligned_fasta/concat/coi_only/COI_AICc/COI_AICc_rooted.treefile") # this is a nexus format. Otherwise, use "read.tree"

# view it to check it
plot(original.coi.p.tree)

# Remove 'parasite' outgroups
no.outgroup.coi.p.tree <- drop.tip(original.coi.p.tree, c('YEWA_941_P_Pr_CTATCCAC-TGTGCGTT', 'YEWA_941_P_Tr_ACGCTTCT-GACTTAGG', 'Proctophyllodes_miliariae'))

# view it to check it
plot(no.outgroup.coi.p.tree)

# Trim tree to one tip per OTU as specified by species delineation (in this case, ABGD)
# first, read in tip labels
no.outgroup.coi.p.tree$tip.label

# drop tips (randomly per OTU)
ABGD.OTU.coi.p.tree <- drop.tip(no.outgroup.coi.p.tree, c('BTNW_765_P_CCAAGTAG-TCACGTTC','YTWA_970_R_GAACCTTC-GTACCTTG','CMWA_952_R_CTAAGACC-CTGTTAGG','CSWA_939_R_TCCATTGC-TCCTACCT','BWWA_412_R_GTACACCT-TTGGTCTC',"BLBW_861_R_AACGCCTT-TTCTCTCG","GWWA_928_R_ACGAGAAC-CGAACTGT","CAWA_366_A_CCGGAATA-CCTACTGA_dd","BTNW_492_R_AAGGACCA-GTGTTCCT","BTBW_850_R_AGGTTCCT-AGGTTCGA","AMRE_967_R_ACGGACTT-TCTCCGAT","CERW_358_A_TCTTCGAC-GAGCAGTA","BPWA_337_C_CAAGCCAA-CTTCGTTC_dd","WEWA_727_R_ACCATCCT-GTTCATGG","PIWA_948_R_ATAACGCC-TTGGTGAG","MOWA_611_R_HC7TNCCX2_subsampled_dd","COYE_868_R_CCAACACT-AGTTCGTC","NAWA_328_A_CCACAACA-CTTGGATG","KEWA_584_R_TCTAGGAG-CGTTGCAA","LOWA_438_A_AACTTGCC-GAAGTTGG"))

# view it
plot(ABGD.OTU.coi.p.tree)

# Change tip labels on this OTU tree
# first, read the remaining tips in to see the correct order needed for relabeling
ABGD.OTU.coi.p.tree$tip.label

# relabel them
new_tiplabels.coi.p.tree <- c("Tyrannidectes_charitomenos", "Tyrannidectes_aff_charitomenos", "Amerodectes_ischyros", "Amerodectes_aff_ischyros", "Amerodectes_helmitheros", "Amerodectes_sp_n_1", "Amerodectes_hribari", "Amerodectes_protonotaria", "Amerodectes_seiurus", "Amerodectes_jonesborensis")

# assign new labels to tree
ABGD.OTU.coi.p.tree$tip.label <- new_tiplabels.coi.p.tree

# plot it to check it
plot(ABGD.OTU.coi.p.tree)

#### Let's export this tree for future use (eg in JANE, empress)
ape::write.tree(ABGD.OTU.coi.p.tree, file = "ABGD.OTU.coi.p.tree.txt")
ape::write.nexus(ABGD.OTU.coi.p.tree, file = "Mites_coi.nex")
ape::write.tree(ABGD.OTU.coi.p.tree, file = "Mites_coi.nwk")



#### ~~ HOST TREE - BirdTree 95 Consensus ####
#### Did not end up using this tree, but including for comparison

# Read in 'host' tree
original.h.tree<- read.nexus("hosts/with_outgroups/consensus95_hosts_ogs.tre") # may also use read.tree if needed

# view it
plot(original.h.tree)

# Remove outgroups
no.outgroup.h.tree <- drop.tip(original.h.tree, c('Icteria_virens', 'Xenoligea_montana', 'Zeledonia_coronata'))

# view it
plot(no.outgroup.h.tree)

# Change tip labels on this host tree (some have old genus names)
# First, get order
no.outgroup.h.tree$tip.label

# change them
new_tiplabels.h.tree <- c('Seiurus_aurocapilla','Helmitheros_vermivorum','Leiothlypis_ruficapilla','Vermivora_cyanoptera','Vermivora_chrysoptera',"Geothlypis_trichas","Geothlypis_formosa","Oporornis_agilis","Geothlypis_philadelphia","Mniotilta_varia","Parkesia_noveboracensis","Parkesia_motacilla","Cardellina_canadensis","Setophaga_citrina","Setophaga_magnolia","Setophaga_cerulea","Setophaga_americana","Setophaga_tigrina","Setophaga_caerulescens","Setophaga_ruticilla","Setophaga_dominica","Setophaga_pinus","Setophaga_discolor","Setophaga_virens","Setophaga_fusca","Setophaga_pensylvanica","Setophaga_striata","Setophaga_castanea","Limnothlypis_swainsonii","Protonotaria_citrea")

# assign new labels
no.outgroup.h.tree$tip.label <- new_tiplabels.h.tree

# view it to check
plot(no.outgroup.h.tree)

#### Some weird stuff here being the CONW placement and the major polytomies



#### ~~ HOST TREE - LOVETTE ET AL. 2010 ####

# Read in tree
original.lovette.h.tree <- read.nexus("/Users/alix22792/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_phylo/Jetstream/hosts/Lovetteetal2010parulids.nex") # may also use read.tree if needed

# view it
plot(original.lovette.h.tree)

# get order of tip labels
original.lovette.h.tree$tip.label

# Remove host outgroups and drop uneeded tips
bird_taxa_to_exclude <- c("'Myioboruscastaneocapilla'","'Myioboruscardonai'","'Myioborusalbifacies'","'Myioboruspariae'","'Myioborustorquatus'","'Myioborusornatus'","'Myioborusmelanocephalus'","'Myioborusalbifrons'","'Myioborusflavivertex'","'Myioborusbrunniceps'","'Myioborusminiatus'","'Myioboruspictus'","'Cardellinaruber'","'Cardellinaversicolor'","'Cardellinarubrifrons'","'Cardellinapusilla'","'Basileuterushypoleucus'","'Basileuterusculicivorus'","'Basileuterustrifasciatus'","'Basileuterustristriatus'","'Basileuterusmelanogenys'","'Basileuterusbelli'","'Basileuterusrufifrons'","'Basileuteruslachrymosa'","'Myiothlypissignatus'","'Myiothlypisnigrocristatus'","'Myiothlypisfulvicauda'","'Myiothlypisrivularis'","'Myiothlypisflaveolus'","'Myiothlypisleucoblepharus'","'Myiothlypisleucophrys'","'Myiothlypisluteoviridis'","'Myiothlypiscoronatus2'","'Myiothlypiscoronatus1'","'Myiothlypisfraseri'","'Myiothlypisconspicillatus'","'Myiothlypiscinereicollis'","'Myiothlypischrysogaster'","'Myiothlypisbivittatus'","'Myiothlypisroraimae'","'Setophagatownsendi'","'Setophagaoccidentalis'","'Setophagachrysoparia'","'Setophagagraciae'","'Setophaganigrescens'","'Setophagasubita'","'Setophagadelicata'","'Setophagaadelaidae'","'Setophagavitellina'","'Setophagapityophila'","'Setophagacoronata'","'Setophagapalmarum'","'Setophagapetechia'","'Setophagapitiayumi'","'Setophagakirtlandii'","'Setophagaangelae'","'Setophagapharetra'","'Setophagaplumbea'","'Setophagabishopi'","'Geothlypisflavovelata'","'Geothlypisnelsoni'","'Geothlypisrostrata'","'Geothlypisbeldingi'","'Geothlypissemiflava'","'Geothlypisspeciosa'","'Geothlypispoliocephala'","'Geothlypisaequinoctialis'","'Geothlypistolmiei'","'Geothlypissemperi'","'Oreothlypisvirginiae'","'Oreothlypisluciae'","'Oreothlypiscrissalis'","'Oreothlypiscelata'","'Oreothlypisperegrina'","'Oreothlypissuperciliosa'","'Oreothlypisgutturalis'","'Vermivorabachmanii'","'Granatelluspelzelni'","'Coerebaflaveola'","'Icteriavirens'","'Microligeapalustris'","'Xenoligeamontana'","'Zeledoniacoronata'","'Spindaliszena'","'Teretistrisfernandinae'")

# assign it
no.outgroups.no.extras.h.tree <- drop.tip(original.lovette.h.tree, bird_taxa_to_exclude)

# plot it to view it
plot(no.outgroups.no.extras.h.tree)


# Change tip labels on the host tree, get order first
no.outgroups.no.extras.h.tree$tip.label

# Create new labels in the correct order
new_tiplabels.no.extras.h.tree <- c("Cardellina_canadensis","Setophaga_virens","Setophaga_discolor","Setophaga_dominica","Setophaga_pinus","Setophaga_caerulescens","Setophaga_pensylvanica","Setophaga_striata","Setophaga_castanea","Setophaga_fusca","Setophaga_americana","Setophaga_cerulea","Setophaga_tigrina","Setophaga_ruticilla","Setophaga_magnolia","Setophaga_citrina","Geothlypis_trichas","Geothlypis_formosa","Geothlypis_philadelphia","Oporornis_agilis","Leiothlypis_ruficapilla","Limnothlypis_swainsonii","Protonotaria_citrea","Mniotilta_varia","Vermivora_chrysoptera","Vermivora_cyanoptera","Parkesia_motacilla","Parkesia_noveboracensis","Helmitheros_vermivorum","Seiurus_aurocapilla")

# Assign new labels
no.outgroups.no.extras.h.tree$tip.label <- new_tiplabels.no.extras.h.tree

# plot it to view it
plot(no.outgroups.no.extras.h.tree)


# Let's export this tree for future use (eg in JANE, empress)
ape::write.tree(no.outgroups.no.extras.h.tree, file = "Lovette_hosts.nwk")
ape::write.nexus(no.outgroups.no.extras.h.tree, file = "Lovette_hosts.nex")


#### ANALYSES: CONCAT AICc (MITES) AND LOVETTE (HOSTS) ####

#### ++ HOST-PARASITE MATRIX ####
#### IMPORTANT: make sure labels match with every tree

# Read in matrix
hpmatrix<-as.matrix(read.table("amero_matrix_newnames.txt", header=T))

# Calculate patristic distances and rotate headings to fit matrix
# 'Parasites' - pick your parasite tree at this point, if you have multiple
parasite.patristic.dist<-as.matrix(cophenetic(ABGD.OTU.p.tree))
parasite.patristic.dist <- parasite.patristic.dist[colnames(hpmatrix), colnames(hpmatrix)]

# Hosts - pick your host tree at this point, if you have multiple
host.patristic.dist<-as.matrix(cophenetic(no.outgroups.no.extras.h.tree))
host.patristic.dist <- host.patristic.dist[rownames(hpmatrix), rownames(hpmatrix)]


#### ++ PARAFIT ####

# Run ParaFit (one run)
# setting the seed helps ensure reproducibility
set.seed(1234)
runParaFit<- function(i){
    tmp <- parafit(host.patristic.dist, parasite.patristic.dist, hpmatrix, nperm=999, test.links=T, correction='cailliez')
    return (tmp)
}

# alternative: correction='lingoes'


# Run parallel ParaFit analyses
cat("Scheduled",length(no.outgroups.no.extras.h.tree),"jobs to run in parallel ... Running ... Nothing will be printed on screen until all runs are completed ...",sep=" ")


# Loop multiple ParaFit Runs (n = 100 in this case)
# Ideal to run multiple runs (> 1) because the results vary slightly over multiple runs given the randomization. Also a good idea to set the seed here so the loop will be repeatable if needing to rerun
res<- lapply(1:100, function(i){
    simplify2array(mclapply(1,runParaFit,mc.preschedule=T,mc.cores=10,mc.set.seed=4321))
})

# View all runs and global p.values (all separate 'objects' from one another
res

# data.class should be 'list'
data.class(res) # SHOULD BE 'LIST'

# Extract p.global values from 'res'
res_p<-sapply(res, "[[",2)

# to see output
res_p

# check class, should be numeric
data.class(res_p)

# make this a matrix for the next step
as.matrix(res_p)

# calculate and view the average p.global value
p_global_mean<-mean(res_p[[1]])
p_global_mean


# Now, need to combine all the parafit runs and p.global values into 1 'object' in order to extract p-values (which then can be used to get the (overall) mean p.global as well as the mean p-values for each link over multiple parafit runs). Also this is useful for making sure the order that is listed in the parafit output is the same as the order the p-values appear (because R likes to reorder these objects in strange ways... Which can lead to the wrong links being associated with 'significance')
res_link<-c()
for (i in 1:100){
    res_link<-cbind(res_link,res[[i]])
}

# to see output
res_link

# should be matrix, all values from 'res' are combined into 'res_link'
data.class(res_link)

# VERY IMPORTANT: Check the order of ParaFit output and p-values
links_parafit_order=cbind(rownames(hpmatrix)[res_link[,1]$link.table[,1]],colnames(hpmatrix)[res_link[,1]$link.table[,2]])

# IMPORTANT: Check if R correctly ordered output (links to p-values)
links_parafit_order

# Create matrix of p values; columns are h-p links, rows are different ParaFit runs
# [certainly, instead of the loop, one can use lapply etc here]
p_table=NULL
for (i in 1:ncol(res_link)){
    p.F1=res_link[,i]$link.table[,4] #col4= p.F1, col6=p.F2
    p_table=rbind(p_table,p.F1)
    #cat(p.F1, "***\n\n")
}
p_table


# Adjust those p-values; columns are hp links, rows are different parafit runs
p_table_adj<-apply(p_table,2,p.adjust,method="BH")

# Check output
p_table_adj

# Should be matrix
data.class(p_table_adj)

# Adjusted p-value means for each individual link – really the only table needed...
p_means<-colMeans(p_table_adj)
p_means_table<-cbind(links_parafit_order,(p_means))
p_means_table_df<-as.data.frame(p_means_table)
p_means_table_df

# be sure to change file name for certain h-p trees used
write.table(p_means_table_df, file='adj_p_each_link_lovette_20220207_setseed1234_newnames.txt',sep='\t', row.names=F, col.names=T,quote=F)


# p-value table, if you want it
ptable<-cbind(links_parafit_order,t(p_table))
colnames(ptable)<-c("Host","Parasite",sprintf("p.F1_tree_%s",seq(1:ncol(res_link))))
ptable

# be sure to change file name for certain h-p trees used
write.table(ptable, file='p_value_table_lovette_20220207_setseed1234_newnames.txt', sep='\t', row.names=F, col.names=T,quote=F)

# p-value talbe with adjusted values per ParaFit run, if you want it
p_tableadj<-cbind(links_parafit_order,t(p_table_adj))
colnames(p_tableadj)<-c("Host","Parasite",sprintf("p.F1.adj_tree_%s",seq(1:ncol(res_link))))

# be sure to change file name for certain h-p trees used
write.table(p_tableadj, file='p_tableadj_lovette_20220207_setseed1234_newnames.txt', sep='\t', row.names=F, col.names=T,quote=F)


#### ++ PLOTTING ####
# Adjust colors and line types
p_colors=NULL

for (i in 1:ncol(p_table_adj)){
    p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
    color="#E69F00" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES
    if(p_min>p_level && p_max>p_level){color="gray75"} # NOT-SIGNIFICANT
    if(p_min<=p_level && p_max<=p_level){color="#0072B2"} # SIGNIFICANT
    p_colors[i]=color
}

p_type=NULL

for (i in 1:ncol(p_table_adj)){
    p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
    lty="solid" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES, can use dotdash
    if(p_min>p_level && p_max>p_level){lty="solid"} # NOT-SIGNIFICANT, can use longdash
    if(p_min<=p_level && p_max<=p_level){lty="solid"} # SIGNIFICANT
    p_type[i]=lty
}

#### Tanglegram to minimize link crossing

ph1 <- no.outgroups.no.extras.h.tree
ph2 <- ABGD.OTU.p.tree
assoc_links <- links_parafit_order
cophy <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

# set as PDF and size of printing area
pdf("ParaFit_only_Lovette_2_20220207_setseed1234_newnames.pdf", height = 30, width = 45)
cophyloplot_min<-plot(cophy, link.lty=p_type, link.col=p_colors, link.lwd=8, fsize=3.5, lwd = 6, pts=F)
dev.off()

# NOTE: If label names are too long.... do this....

# EMMANUEL PARADIS AT http://grokbase.com/t/r/r-sig-phylo/12c3djkx3h/modification-of-figures-produced-using-cophyloplot
# FIRST MAKE A COPY OF THE INTERNAL FUNCTION
# titi<-plotCophylo2
# then do fix(titi) and modify what you want (here look for calls to
# text(... cex = ...), save and close. Do fix(cophyloplot) and change
# "plotCophylo2(..." by "titi(...", save and close.
# fix(titi) #
# if (show.tip.label) {
#   text(a[1:N.tip.x, ], cex = 0.5, font = font, pos = 4, labels = x$tip.label)
#   text(b2[1:N.tip.y, ], cex = 0.5, font = font, pos = 2,
#        labels = y$tip.label)
# }
#
# titi
# fix(cophyloplot) # Do fix(cophyloplot) and change "plotCophylo2(..." by "titi(...", save and close.






#### ++ PACo - 'old' way ####
#### useful for calculating m2xy contributions on jackknife

# PACo function code.. The PACo function, defined below, transforms the host and parasite distance matrices into the respective matrices of Principal Coordinates (pcoa of ape) and duplicates taxa (if necessary) to accommodate multiple host-parasite associations:
PACo_old <- function (H.dist, P.dist, HP.bin)
{
    HP.bin <- which(HP.bin > 0, arr.in=TRUE)
    H.PCo <- pcoa(H.dist, correction="cailliez")$vectors #Performs PCo of Host distances
    P.PCo <- pcoa(P.dist, correction="cailliez")$vectors #Performs PCo of Parasite distances
    H.PCo <- H.PCo[HP.bin[,1],] #adjust Host PCo vectors
    P.PCo <- P.PCo[HP.bin[,2],]  ##adjust Parasite PCo vectors
    list (H.PCo = H.PCo, P.PCo = P.PCo)
}


# Calculate patristic distances
# Pick your h-p trees as this point
mites.paco.D <- cophenetic(ABGD.OTU.p.tree) # symbionts
hosts.paco.D <- cophenetic(no.outgroups.no.extras.h.tree) # hosts

# rescale hpmatrix to cophenetic trees
mites.paco.D <- mites.paco.D[colnames(hpmatrix), colnames(hpmatrix)]
hosts.paco.D <- hosts.paco.D[rownames(hpmatrix), rownames(hpmatrix)]

# apply PACo to the data
PACo.fit <- PACo_old(hosts.paco.D, mites.paco.D, hpmatrix)
PACo.procru <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)
NLinks = sum(hpmatrix)

# a warning may be produced which indicates that the host input matrix has fewer columns than the parasite counterpart. No action by the user is required since the narrower matrix is completed with zero columns


# Ordination plot stuff, if desired
# host ordination matrix
HostX <- PACo.procru$X

# Parasite ordination matrix, scaled and rotated to fit HostX
ParY <- PACo.procru$Yrot

# Plotting host and parasite ordinations
plot(HostX, asp=1, pch=46)
points(ParY, pch=1)
arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15, xpd=FALSE)
HostX <- unique(PACo.procru$X)

# unique() removes duplicated points - convenient for labelling of points below
ParY <- unique(PACo.procru$Yrot)

# interactive labelling
identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8)
identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)




# Goodness of fit: the following code computes the residual sum of squares m2xy and performs a randomization of the host-parasite association matrix to establish the probability (P) under H0

# observed sum of squares
m2.obs <- PACo.procru$ss

# number of permutations, this can be as low as 10,000 to 100,000
N.perm = 100000

# p-value
P.value = 0

# set a seed either randomly or reproducibly
# seed <-.Random.seed[trunc(runif(1,1,626))]
# set.seed(seed)
set.seed(1234) ### use this option to obtain reproducible randomizations



# IMPORTANT: !!! Given that append=TRUE (below), the file created (m2_perm.txt) should be deleted or renamed prior to a new analysis. Otherwise the values generated in the new run will be appended to those produced in the previous one.

for (n in c(1:N.perm))
{
    if (NLinks <= nrow(hpmatrix) | NLinks <= ncol(hpmatrix)) 	#control statement to avoid all parasites beig associated to a single host
    {	flag2 <- TRUE
    while (flag2 == TRUE)	{
        HP.perm <- t(apply(hpmatrix,1,sample))
        if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
    }
    } else { HP.perm <- t(apply(hpmatrix,1,sample))} #permutes each HP row independently
    PACo.perm <- PACo_old(hosts.paco.D, mites.paco.D, HP.perm)
    m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
    write (m2.perm, file = "PACo_m2perm_concatmites_lovette_setsed1234_newnames_xxx.txt", sep ="\t", append =TRUE) #option to save m2 from each permutation
    if (m2.perm <= m2.obs)
    {P.value = P.value + 1}
}


# To conclude the goodness-of-fit test:
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P value = ", P.value, " based on ", N.perm,"permutations.")

# warnings are okay, no effect on analysis
# so m2 = whatever, and P-value = something. IF... P = 10^-5, then in only 1 out of 100,000 permutations (1/100,000 = 10^-5), the residual sum of squares was smaller than "m2", and congruence between the host and parasite phylogenies is statistically significant at the a = 0.05 level.


# Jackknife procedure
# this will produce a bar chart of squared residuals... the links below the abline (which should be a horizontal dotted line/median global squared residual) mean they contribute relatively little to the m2xy and thus likely represent coevolutionary links. the larger the squared residuals, the less they are contributing to the cophylogenetic pattern (adding noise or showing discordant phylogenetic patterns)... but should consider their confidence intervals (if the CIs are quite broad, then it's difficult to evaluate their real contribution to the cophylogenetic pattern observed)

HP.ones <- which(hpmatrix > 0, arr.in=TRUE)

# empty matrix of jackknifed squared residuals
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)

# colnames identify the H-P link
colnames (SQres.jackn) <- paste(rownames(PACo.procru$X),rownames(PACo.procru$Yrot), sep="-")

# needed to compute 95% confidence intervals.
t.critical = qt(0.975,NLinks-1)

# PACo setting the ith link = 0
for(i in c(1:NLinks))
{HP.ind <- hpmatrix
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo_old(hosts.paco.D, mites.paco.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)

# Append residuals to matrix of jackknifed squared residuals
SQres.jackn [i, ] <- res.Proc.ind
}

# Jackknifed residuals are squared
SQres.jackn <- SQres.jackn**2

# Vector of original square residuals
SQres <- (residuals (PACo.procru))**2


# Jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) # apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) # mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) # standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) # upper 95% confidence interval

# plotting
pdf("PACo_jackknife_concatmites_lovette_setset1234_newnames.pdf", width = 14, height = 14)

mar.default <- c(4,4,1,2) + 0.1 # adjusting margins to help with the long labels
par(mar = mar.default + c(17, 1, 1, 5)) # continue to adjust margins to help with long labels

pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.5, col="gray", xlab= "", ylab= "Squared residuals", ylim=c(0, max(phi.UCI)), cex.lab=2)
mtext("Host-mite links", cex=2, side = 1, line = 17) # line = should match the first value in mar.default. If you put this label in xlab above, it will be placed in a strange spot
text(pat.bar, par("usr")[3] - 0.001, srt = 300, adj = -0.00005, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.75)  # can adjust all of these values to change how the labels are displayed
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2, xpd=FALSE)

dev.off()


#### ++ PACo - 'new way' ####

# Transform phylogenies to distance matrices.
# Pick your trees as this point
mites.paco.D <- cophenetic(ABGD.OTU.p.tree) # symbionts
hosts.paco.D <- cophenetic(no.outgroups.no.extras.h.tree) # hosts

# Bundle the data together and add Principal Coordinates
D <- prepare_paco_data(H = hosts.paco.D, P = mites.paco.D, HP = hpmatrix)
D <- add_pcoord(D, correction = 'cailliez')

# Perform cophylogenetic analyses
D <-  paco::PACo(D, nperm = 999, seed = 12, method = "r0", symmetric = FALSE, proc.warnings = FALSE, shuffled = FALSE)

# symmetric = FALSE assumes that the parasite tracks the evolution of the host. If set to TRUE, this means that both phylogenies are standardized prior to superimposition resulting in the best-fit of the superimposition is independent of both phylos (one does not track the other). The method = 'r0' has to do with how the association matrix is permuted (the randomisation algorithms). r0 maintains the row sums as well as the # of interactions per species and should be used under the assumption that the column group (in my matrix, the mites) tracks the evolution of the row (hosts) group, therefore maintaining the degree of the species in the row group. Could try a few methods here.

print(D$gof) # goodness of fit, showing p-value, m2xy value, and n

# return interaction-specific cophylogenetic contributions based on a jackknifing procedure
D <- paco::paco_links(D)

print(D$jackknife) # bias-corrected residuals for each link.

# return interaction-specific cophylogenetic contributions in terms of individual superimposition residuals
paco_resids <- paco::residuals_paco(D$proc)

print(paco_resids) # Procrustes residuals are returned. Smaller residual distance = stronger cophylo signal

# to visualise this we use the ape function cophyloplot weighted by interaction contribution
# first we must make a list out of the interaction matrix
assoc <- data.frame(host=rownames(hpmatrix)[which(hpmatrix==1, arr.ind=TRUE)[,'row']], mite=colnames(hpmatrix)[which(hpmatrix==1, arr.ind=TRUE)[,'col']])

# to weight the interactions we use the cophylogenetic contribution transformed to best show the differences graphically (inversely proportional to the resids = thicker lines = interaction that shows stronger support for phylogenetic congruence)
weight <- (paco_resids^-2)/50

cophyloplot(no.outgroups.no.extras.h.tree, ABGD.OTU.p.tree, assoc, show.tip.label=FALSE, use.edge.length=FALSE, lwd=weight, col='steelblue', length.line=0, gap=-20, space=60)

cophyloplot(no.outgroups.no.extras.h.tree, ABGD.OTU.p.tree, assoc, show.tip.label=T, use.edge.length=FALSE, lwd=weight, col='steelblue', length.line=0, gap=20, space=100, rotate=F)


# Let's combine the weights of the links from PACo to the tanglegram from ParaFit
# check parafit order
links_parafit_order

# check paco weights order
weight

# they are not in the same order. So make the weights ordered as ParaFit (because colors are ordered by ParaFit and minimize the link crossing)

paco_weight_parafit_order_man <- c(2.677715, 3.379930, 4.913892, 3.331220, 3.061097, 5.565896, 4.320257, 4.424552, 7.922728, 8.598833, 7.047518, 6.459642, 3.257489, 5.120214, 6.458951, 6.698400, 5.894914, 5.544560, 6.142252, 5.814481, 7.394063, 6.968635, 3.327400, 10.724615, 8.581129, 5.027256, 7.079873, 7.701262, 12.795780, 13.762410, 3.639649, 5.701532)

ph1 <- no.outgroups.no.extras.h.tree
ph2 <- ABGD.OTU.p.tree
assoc_links <- links_parafit_order
cophy2 <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

plot(cophy2)

# Figure in publication
pdf("PACo_ParaFit_combo_20220505_Lovette.pdf", height = 30, width = 45) # sets size of printing area
cophyloplot_min<-plot(cophy2, link.lty=p_type, link.col=p_colors, link.lwd=(paco_weight_parafit_order_man)*2, fsize=5,pts=F, lwd=10)
dev.off()









#### ANALYSES: COI AICc (MITES) AND LOVETTE (HOSTS) ####

#### ++ HOST-PARASITE MATRIX ####
#### IMPORTANT: make sure labels match with every tree

# Read in matrix
hpmatrix<-as.matrix(read.table("amero_matrix_newnames.txt", header=T))

# Calculate patristic distances and rotate headings to fit matrix
# 'Parasites' - pick your parasite tree at this point, if you have multiple
parasite.patristic.dist<-as.matrix(cophenetic(ABGD.OTU.coi.p.tree))
parasite.patristic.dist <- parasite.patristic.dist[colnames(hpmatrix), colnames(hpmatrix)]


# Hosts - pick your host tree at this point, if you have multiple
host.patristic.dist<-as.matrix(cophenetic(no.outgroups.no.extras.h.tree))
host.patristic.dist <- host.patristic.dist[rownames(hpmatrix), rownames(hpmatrix)]





#### ++ PARAFIT  ####
# Run ParaFit (one run)
# setting the seed helps ensure reproducibility
set.seed(1234)
runParaFit<- function(i){
    tmp <- parafit(host.patristic.dist, parasite.patristic.dist, hpmatrix, nperm=999, test.links=T, correction='cailliez')
    return (tmp)
}

# alternative: correction='lingoes'


# Run parallel ParaFit analyses
cat("Scheduled",length(no.outgroups.no.extras.h.tree),"jobs to run in parallel ... Running ... Nothing will be printed on screen until all runs are completed ...",sep=" ")


# Loop multiple ParaFit Runs (n = 100 in this case)
# Ideal to run multiple runs (> 1) because the results vary slightly over multiple runs given the randomization. Also a good idea to set the seed here so the loop will be repeatable if needing to rerun
res<- lapply(1:100, function(i){
    simplify2array(mclapply(1,runParaFit,mc.preschedule=T,mc.cores=10,mc.set.seed=F))
})

# View all runs and global p.values (all separate 'objects' from one another)
res

# data.class should be 'list'
data.class(res)


# Extract p.global values from 'res'
res_p<-sapply(res, "[[",2)

# to see output
res_p

# check class, should be numeric
data.class(res_p)

# make this a matrix for the next step
as.matrix(res_p)

# calculate and view the average p.global value
p_global_mean<-mean(res_p[[1]])
p_global_mean


# Now, need to combine all the parafit runs and p.global values into 1 'object' in order to extract p-values (which then can be used to get the (overall) mean p.global as well as the mean p-values for each link over multiple parafit runs). Also this is useful for making sure the order that is listed in the parafit output is the same as the order the p-values appear (because R likes to reorder these objects in strange ways... Which can lead to the wrong links being associated with 'significance')
res_link<-c()
for (i in 1:100){
    res_link<-cbind(res_link,res[[i]])
}

# to see output
res_link

# should be matrix, all values from 'res' are combined into 'res_link'
data.class(res_link)

# VERY IMPORTANT: Check the order of ParaFit output and p-values
links_parafit_order=cbind(rownames(hpmatrix)[res_link[,1]$link.table[,1]],colnames(hpmatrix)[res_link[,1]$link.table[,2]])

# IMPORTANT: Check if R correctly ordered output (links to p-values)
links_parafit_order


# Create matrix of p values; columns are h-p links, rows are different ParaFit runs
# [certainly, instead of the loop, one can use lapply etc here]
p_table=NULL
for (i in 1:ncol(res_link)){
    p.F1=res_link[,i]$link.table[,4] #col4= p.F1, col6=p.F2
    p_table=rbind(p_table,p.F1)
    #cat(p.F1, "***\n\n")
}
p_table


# Adjust those p-values; columns are hp links, rows are different ParaFit runs
p_table_adj<-apply(p_table,2,p.adjust,method="BH")

# Check output
p_table_adj

# Should be matrix
data.class(p_table_adj)


# Adjusted p-value means for each individual link – really the only table needed...
p_means<-colMeans(p_table_adj)
p_means_table<-cbind(links_parafit_order,(p_means))
p_means_table_df<-as.data.frame(p_means_table)
p_means_table_df

# be sure to change file name for certain h-p trees used
write.table(p_means_table_df, file='adj_p_each_link_lovette_COI_20220505_setseed1234.txt',sep='\t', row.names=F, col.names=T,quote=F)


# p-value table, if you want it
ptable<-cbind(links_parafit_order,t(p_table))
colnames(ptable)<-c("Host","Parasite",sprintf("p.F1_tree_%s",seq(1:ncol(res_link))))
ptable

# be sure to change file name for certain h-p trees used
write.table(ptable, file='p_value_table_lovette_COI_20220505_setseed1234.txt', sep='\t', row.names=F, col.names=T,quote=F)


# p-value table with adjusted values per ParaFit run, if you want it
p_tableadj<-cbind(links_parafit_order,t(p_table_adj))
colnames(p_tableadj)<-c("Host","Parasite",sprintf("p.F1.adj_tree_%s",seq(1:ncol(res_link))))

# be sure to change file name for certain h-p trees used
write.table(p_tableadj, file='p_tableadj_lovette_COI_20220505_setseed1234.txt', sep='\t', row.names=F, col.names=T,quote=F)


#### ++ PLOTTING ####
# Adjust colors and line types
p_colors=NULL

for (i in 1:ncol(p_table_adj)){
    p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
    color="#E69F00" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES
    if(p_min>p_level && p_max>p_level){color="gray75"} # NOT-SIGNIFICANT
    if(p_min<=p_level && p_max<=p_level){color="#0072B2"} # SIGNIFICANT
    p_colors[i]=color
}

p_type=NULL

for (i in 1:ncol(p_table_adj)){
    p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
    lty="solid" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES
    if(p_min>p_level && p_max>p_level){lty="solid"} # NOT-SIGNIFICANT
    if(p_min<=p_level && p_max<=p_level){lty="solid"} # SIGNIFICANT
    p_type[i]=lty
}

#### Tanglegram to minimize link crossing
ph1 <- no.outgroups.no.extras.h.tree
ph2 <- ABGD.OTU.coi.p.tree
assoc_links <- links_parafit_order
cophy <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

pdf("ParaFit_only_Lovette_COI_AICc_20220505_setseed1234_newnames.pdf", height = 30, width = 45) # sets size of printing area
cophyloplot_min<-plot(cophy, link.lty=p_type, link.col=p_colors, link.lwd=8, fsize=3.5,pts=F)
dev.off()


#### ++ PACo - 'old' way ####
#### useful for calculating m2xy contributions on jackknife

# PACo function code.. The PACo function, defined below, transforms the host and parasite distance matrices into the respective matrices of Principal Coordinates (pcoa of ape) and duplicates taxa (if necessary) to accommodate multiple host-parasite associations:
PACo_old <- function (H.dist, P.dist, HP.bin)
{
    HP.bin <- which(HP.bin > 0, arr.in=TRUE)
    H.PCo <- pcoa(H.dist, correction="cailliez")$vectors #Performs PCo of Host distances
    P.PCo <- pcoa(P.dist, correction="cailliez")$vectors #Performs PCo of Parasite distances
    H.PCo <- H.PCo[HP.bin[,1],] #adjust Host PCo vectors
    P.PCo <- P.PCo[HP.bin[,2],]  ##adjust Parasite PCo vectors
    list (H.PCo = H.PCo, P.PCo = P.PCo)
}

# Calculate patristic distances
# Pick your h-p trees as this point
mites.paco.D <- cophenetic(ABGD.OTU.coi.p.tree) # symbionts
hosts.paco.D <- cophenetic(no.outgroups.no.extras.h.tree) # hosts

#### rescale hpmatrix to cophenetic trees
mites.paco.D <- mites.paco.D[colnames(hpmatrix), colnames(hpmatrix)]
hosts.paco.D <- hosts.paco.D[rownames(hpmatrix), rownames(hpmatrix)]

# apply PACo to the data
PACo.fit <- PACo_old(hosts.paco.D, mites.paco.D, hpmatrix)
PACo.procru <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)
NLinks = sum(hpmatrix)

# a warning will be produced which indicates that the host input matrix has fewer columns than the parasite counterpart. No action by the user is required since the narrower matrix is completed with zero columns


# Ordination plot stuff, if desired
# host ordination matrix
HostX <- PACo.procru$X
# parasite ordination matrix, scaled and rotated to fit HostX
ParY <- PACo.procru$Yrot

# Plotting host and parasite ordinations
plot(HostX, asp=1, pch=46)
points(ParY, pch=1)
arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15, xpd=FALSE)
HostX <- unique(PACo.procru$X)

# unique() removes duplicated points - convenient for labelling of points below
ParY <- unique(PACo.procru$Yrot)

# interactive labelling
identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8)
identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)




# Goodness of fit: the following code computes the residual sum of squares m2xy and performs a randomization of the host-parasite association matrix to establish the probability (P) under H0

# observed sum of squares
m2.obs <- PACo.procru$ss

# this can be as low as 10,000 to 100,000
N.perm = 100000
P.value = 0

# set a seed either randomly or reproducibly
# seed <-.Random.seed[trunc(runif(1,1,626))]
# set.seed(seed)
set.seed(1234) ### use this option to obtain reproducible randomizations



# IMPORTANT: !!! Given that append=TRUE (below), the file created (m2_perm.txt) should be deleted or renamed prior to a new analysis. Otherwise the values generated in the new run will be appended to those produced in the previous one.

for (n in c(1:N.perm))
{
    if (NLinks <= nrow(hpmatrix) | NLinks <= ncol(hpmatrix)) 	#control statement to avoid all parasites beig associated to a single host
    {	flag2 <- TRUE
    while (flag2 == TRUE)	{
        HP.perm <- t(apply(hpmatrix,1,sample))
        if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
    }
    } else { HP.perm <- t(apply(hpmatrix,1,sample))} #permutes each HP row independently
    PACo.perm <- PACo_old(hosts.paco.D, mites.paco.D, HP.perm)
    m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
    write (m2.perm, file = "PACo_m2perm_COImites_lovette_xxx.txt", sep ="\t", append =TRUE) #option to save m2 from each permutation
    if (m2.perm <= m2.obs)
    {P.value = P.value + 1}
}



# To conclude the goodness-of-fit test:
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P value = ", P.value, " based on ", N.perm,"permutations.")

# warnings are okay, no effect on analysis
# so m2 = whatever, and P-value = something. IF... P = 10^-5, then in only 1 out of 100,000 permutations (1/100,000 = 10^-5), the residual sum of squares was smaller than "m2", and congruence between the host and parasite phylogenies is statistically significant at the a = 0.05 level.


# Jackknife procedure
# this code will produce a bar chart of squared residuals... the links below the abline (which should be a horizontal dotted line/median global squared residual) mean they contribute relatively little to the m2xy and thus likely represent coevolutionary links. the larger the squared residuals, the less they are contributing to the cophylogenetic pattern (adding noise or showing discordant phylogenetic patterns)... but should consider their confidence intervals (if the CIs are quite broad, then it's difficult to evaluate their real contribution to the cophylogenetic pattern observed)

HP.ones <- which(hpmatrix > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(PACo.procru$X),rownames(PACo.procru$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- hpmatrix
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo_old(hosts.paco.D, mites.paco.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
}

SQres.jackn <- SQres.jackn**2 #Jackknifed residuals are squared
SQres <- (residuals (PACo.procru))**2 # Vector of original square residuals


#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) #standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) #upper 95% confidence interval


# plotting
pdf("PACo_jackknife_COImites_lovette.pdf", width = 14, height = 14)

mar.default <- c(4,4,1,2) + 0.1 # adjusting margins to help with the long labels
par(mar = mar.default + c(17, 1, 1, 5)) # continue to adjust margins to help with long labels

pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.5, col="gray", xlab= "", ylab= "Squared residuals", ylim=c(0, max(phi.UCI)), cex.lab=2)
mtext("Host-mite links", cex=2, side = 1, line = 17) # line = should match the first value in mar.default. If you put this label in xlab above, it will be placed in a strange spot
text(pat.bar, par("usr")[3] - 0.001, srt = 300, adj = -0.00005, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.75)  # can adjust all of these values to change how the labels are displayed
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2, xpd=FALSE)

dev.off()





#### ++ PACO - 'new way' ####

#### Transform phylogenies to distance matrices.
# Pick your trees as this point
mites.paco.D <- cophenetic(ABGD.OTU.coi.p.tree) # symbionts
hosts.paco.D <- cophenetic(no.outgroups.no.extras.h.tree) # hosts

# Bundle the data together and add Principal Coordinates
D <- prepare_paco_data(H = hosts.paco.D, P = mites.paco.D, HP = hpmatrix)
D <- add_pcoord(D, correction = 'cailliez')

# Perform cophylogenetic analyses
D <-  paco::PACo(D, nperm = 999, seed = 12, method = "r0", symmetric = FALSE, proc.warnings = FALSE, shuffled = FALSE)

# symmetric = FALSE assumes that the parasite tracks the evolution of the host. If set to TRUE, this means that both phylogenies are standardized prior to superimposition resulting in the best-fit of the superimposition is independent of both phylos (one does not track the other). The method = 'r0' has to do with how the association matrix is permuted (the randomisation algorithms). r0 maintains the row sums as well as the # of interactions per species and should be used under the assumption that the column group (in my matrix, the mites) tracks the evolution of the row (hosts) group, therefore maintaining the degree of the species in the row group. Could try a few methods here.

print(D$gof) # goodness of fit, showing p-value, m2xy value, and n

# return interaction-specific cophylogenetic contributions based on a jackknifing procedure
D <- paco::paco_links(D)

print(D$jackknife) # bias-corrected residuals for each link.

# return interaction-specific cophylogenetic contributions in terms of individual superimposition residuals
paco_resids <- paco::residuals_paco(D$proc)

print(paco_resids) # Procrustes residuals are returned. Smaller residual distance = stronger cophylo signal

# to visualise this we use the ape function cophyloplot weighted by interaction contribution
# first we must make a list out of the interaction matrix
assoc <- data.frame(host=rownames(hpmatrix)[which(hpmatrix==1, arr.ind=TRUE)[,'row']], mite=colnames(hpmatrix)[which(hpmatrix==1, arr.ind=TRUE)[,'col']])

# to weight the interactions we use the cophylogenetic contribution transformed to best show the differences graphically (inversely proportional to the resids = thicker lines = interaction that shows stronger support for phylogenetic congruence)
weight <- (paco_resids^-2)/50

cophyloplot(no.outgroups.no.extras.h.tree, ABGD.OTU.coi.p.tree, assoc, show.tip.label=FALSE, use.edge.length=FALSE, lwd=weight, col='steelblue', length.line=0, gap=-20, space=60)


cophyloplot(no.outgroups.no.extras.h.tree, ABGD.OTU.coi.p.tree, assoc, show.tip.label=T, use.edge.length=FALSE, lwd=weight, col='steelblue', length.line=0, gap=20, space=100, rotate=F)


# Let's combine the weights of the links from PACo to the tanglegram from ParaFit
# check parafit order
links_parafit_order

# check paco weights order
weight

# they are not in the same order. So make the weights ordered as ParaFit (because colors are ordered by ParaFit and minimize the link crossing)

paco_weight_parafit_order_man <- c(3.188860, 2.923011, 5.095384, 3.663939, 3.330818, 4.688039, 3.515545, 3.942642, 6.320274, 6.013151, 10.669670, 9.742959, 3.447013, 5.717407, 6.086983, 6.046664, 5.289228, 6.264536, 5.482961, 5.223819, 6.594282, 7.733841, 3.940018, 14.840464, 9.513190, 4.717104, 6.080592, 6.501371, 8.453783, 9.476344, 3.167626, 4.342135)

ph1 <- no.outgroups.no.extras.h.tree
ph2 <- ABGD.OTU.coi.p.tree
assoc_links <- links_parafit_order
cophy2 <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

plot(cophy2)

# Figure in publication
pdf("PACo_ParaFit_combo_20220505_Lovette_COI_AICc.pdf", height = 30, width = 45) # sets size of printing area
cophyloplot_min<-plot(cophy2, link.lty=p_type, link.col=p_colors, link.lwd=(paco_weight_parafit_order_man)*2, fsize=5,pts=F,lwd=10)
dev.off()


