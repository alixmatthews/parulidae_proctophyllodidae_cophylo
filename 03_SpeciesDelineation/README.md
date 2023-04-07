## Files in this directory:
- `GeneticDistances.R`

Species delination attempts from alignments:

### R

```GeneticDistances.R``` - this calculates genetic distances of species based on .fasta files using multiple algorithms



### ABGD

Installed it on XSEDE: 
```
	$ wget https://bioinfo.mnhn.fr/abi/public/abgd/last.tgz
	$ tar -xvzf last.tgz
	$ cd Abgd 
	$ make
```

Ready to roll!


```
#pwd: /vol_c/mt_genes/20211115/aligned_fasta/MAFFT
/home/alixmatthews/Abgd/abgd -a -d 0 -p 0.0001 -P 0.1 coi_cat.fasta.mafftaligned.fasta -o ABGD

# -d 0 refers to K2P distances
# results are the same as with HOWA. Looks like we have 10 species groups. *charitomenos* and *ischyros* split into 2 groups 

```

```
/home/alixmatthews/Abgd/abgd -a -d 1 -p 0.0001 -P 0.1 coi_cat.fasta.mafftaligned.fasta -o ABGD_JC

# -d 1 refers to Jukes-Cantor distances
# Same results as K2P
```

```
/home/alixmatthews/Abgd/abgd -a -d 3 -p 0.0001 -P 0.1 coi_cat.fasta.mafftaligned.fasta -o ABGD_raw

# -d 3 refers to raw genetic distances
# Same results as K2P and JC
```

### PTP/bPTP

Citations:

**PTP**:
A General Species Delimitation Method with Applications to Phylogenetic Placements. Zhang, Jiajie, Kapli, P., Pavlidis, P., and Stamatakis, A. Bioinformatics (Oxford, England)(2013), 29 (22): 2869-2876

**bPTP**:
A General Species Delimitation Method with Applications to Phylogenetic Placements. Zhang, Jiajie, Kapli, P., Pavlidis, P., and Stamatakis, A. Bioinformatics (Oxford, England)(2013), 29 (22): 2869-2876

Did on the online portal:
IQ-Tree AICc .treefile concatenated

Job ID was 93112; said it was unrooted, asked to remove outgroups. 100,000 MCMC generations, 0.1 burnin (will check for convergence), seed 123, 100 thinning







