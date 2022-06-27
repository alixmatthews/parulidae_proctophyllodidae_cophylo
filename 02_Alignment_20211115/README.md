02_Alignment_20211115 - final and best attempt (all 11 mt genes that were assembled) with the removal of the HOWA sample.

Alright, alright. It is clear the HOWA (*S. citrina*) sample is causing an issue, given the concatenated tree was poorly supported and the topology was really weird, with *charitomenos* and *ischryos* being sister to one another, while the coalescent was looking more appropriate, AND the gene trees showed that HOWA was falling into different clades depending on the particular gene. I also looked at distance gene trees with this sample included to check for numts... if it was a numt, this sample would be very distant from the others (it wasn't). SO, let's remove this HOWA sample from the pile and take a look again.

### Alignment with MAFFT

Alignment by each gene separately:

```
for file in *.fasta;do 
mafft --auto $file >$file.mafftaligned.fasta; 
done
```

- These files are found in the `mafft_fastas` directory here.

Concatenate these individual gene alignments by sample ID:

```
fasconcat-g -l -s
```

- This concatenated file is found in the `mafft_fastas` directory here.

Check it out in SeaView, use amino acid translation (genetic code = 5, invertebrate mitochondrial)

Looks good!


### IQ-TREE

Partition file is each of the genes (i.e., partitioned by gene)

Merit: BIC
```
iqtree -s Mites_concat.fas -spp Mites_concat_partition.txt -m MFP+MERGE -bb 1000 -pre Mites_tree_concat_BIC
mkdir iqtree_BIC
mv Mites_tree_concat.* iqtree_BIC
```
open the .treefile
THIS concatenated tree … actually looks really good!!!


Merit: AIC
```
iqtree -s Mites_concat.fas -spp Mites_concat_partition.txt -m MFP+MERGE -merit AIC -bb 1000 -pre Mites_tree_concat_AIC
mkdir iqtree_AIC
mv Mites_tree_concat_AIC.* iqtree_AIC
```
open the .treefile
This tree looks really good too, higher support than BIC


Merit: AICc
```
iqtree -s Mites_concat.fas -spp Mites_concat_partition.txt -m MFP+MERGE -merit AICc -bb 1000 -pre Mites_tree_concat_AICc
mkdir iqtree_AICc
mv Mites_tree_concat_AICc.* iqtree_AICc
```
open the .treefile
This tree looks really good too, higher support than BIC and AIC. Highest support! << this (AICc) was the final concatenated tree that I moved forward with.



### ASTRAL: went ahead to do gene trees -> coalescent tree again for old time's sake

IQ-TREE by gene
```
for gene in *.fasta; do
iqtree -s $gene -m MFP -bb 1000;
done
```

Did some data rearrangements

```
Mkdir concat_tree
Mv Mites_concat* concat_tree
Mv Mites_tree_concat* concat_tree
Mkdir gene_trees
Mv *.fasta.mafftaligned.fasta.* gene_trees
Cd gene_trees
cat *treefile > mt_genes_trees.trees
```

And then ASTRAL on these gene_trees.trees files

```
astral -i mt_genes_trees.trees -o mite_species_tree.tre
```

This is the final coalescent tree I moved forward with


### Now to compare with COI only...

Tried BIC, AIC, and AICc in the same way, AICc had the highest support, so moved forward with the AICc tree
```
iqtree -s coi_cat.fasta.mafftaligned.fasta -m MFP -merit AICc -bb 1000 -pre COI_AICc
```

Best-fit model according to AICc: TIM3+F+I+G4 ** this tree has the highest support (even though it is still very low for some clades)...



