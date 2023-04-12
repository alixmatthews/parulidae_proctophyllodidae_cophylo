# Host tree and cophylogenetic tests

---
## Files in this directory:
Various input/output files for cophylogenetic analysis 

- `COI_AICc_rooted.treefile`
- `Lovette_Mites_coi.nex`
- `Lovette_Mites_coi_jane_renamed.nex`
- `Lovette_Mites_concat_AICc.nex`
- `Lovette_Mites_concatenated_jane_renamed.nex`
- `Lovetteetal2010parulids.nex`
- `Mites_tree_concat_AICc_rooted.treefile`
- `ParaFit_PACo_Clean.R`
- `amero_matrix.csv`
- `amero_matrix_newnames.txt`
- `consensus95_hosts_ogs.tre`
- `make_tangle.py`

---

## Obtain a host tree

### Option 1. BirdTree Online Database
- Download bird trees by the following: Birdtree.org -> phylogeny subsets -> copy hosts over because some names are different -> email address and source of Trees should be Hackett All Species (a set of 10k trees with 9993 OTUS each) backbone -> 1000 trees created
- Use SumTrees version 4.5.2 to summarize these 1000 trees
- Need to install DendroPy version 4.5.2 first, on XSEDE: 
``` 
sudo pip install -U dendropy 
```


- Tried a maximum-clade credibility tree first (but discarded it in the end because of really low support)

```
sumtrees.py --summary-target=mcct --burnin=200 --support-as-labels --output-tree-filepath=MCC_hosts_ogs.tre parulid_hosts_ogs.nex
```

- Then tried to do a 50% majority-rule consensus tree

```
sumtrees.py --summary-target=consensus --min-clade-freq=0.5 --burnin=200 --support-as-labels --output-tree-filepath=consensus50_hosts_ogs.tre parulid_hosts_ogs.nex
```

- Then tried to do a 95% majority-rule consensus tree (which is more conservative)

```
sumtrees.py --summary-target=consensus --min-clade-freq=0.95 --burnin=200 --support-as-labels --output-tree-filepath=consensus95_hosts_ogs.tre parulid_hosts_ogs.nex
```


### Option 2. Obtain directly from Open Tree of Life and Lovette et al. 2010
- Obtain tree here: `https://tree.opentreeoflife.org/curator/study/view/pg_2591/?tab=home`
- This is the one I ended up using for downstream processess... (`Lovette_Mites_coi.nex`)
- BirdTree tree was not as clear and had some odd sister clades which conflicted with Lovette et al. 2010 (which is what is used for taxonomy today)... 


---

## Cophylogenetic Tests

### ParaFit and PACo

Necessary files:

- `amero_matrix_newnames.txt`: association matrix of mites (columns) and hosts (rows). Make sure names match identically with phylogeny files in order for R script below to run properly. This is the one with the new names based on morphological analyses

- `ParaFit_PACo.R` - R script that makes .nwk files of phylogenies using ape (for import into JANE), also run through the entirety of the ParaFit and PACo analyses (PACo in both the 'old' way, which makes the plots showing the mxy2 value for each association and also the 95% CIs for each link, as well as the 'new and improved' way which is much faster, but (as far as I know) doesn't make those plots)

- `Lovetteetal2010parulids.nex`: Lovette et al. 2010 parulid tree file from the Open Tree of Life

- `consensus95_hosts_ogs.tre`: BirdTree 95% consensus tree (see above section for how this was generated)

- `COI_AICc_rooted.treefile`: mite tree (COI only), rooted, using AICc 

- `Mites_tree_concat_AICc_rooted.treefile`: mite tree (all mt genes, concatenated), rooted, using AICc



### JANE
- `amero_matrix.csv`: association_matrix.txt (old names) converted from .txt to .csv in Excel, necessary to do so for the `make_tangle.py` script below

- `make_tangle.py`: script created by Andrew Sweet to take in .nwk files of phylogenies (can be created in the .R script or however else you want to do it), an association matrix (above) and will create a .nex file necessary for JANE.

Then I ran the following to make the .nex files necessary for import into JANE (unless you want to do it by hand, ha!):

```
python3 make_tangle.py -m amero_matrix.csv -s Lovette_hosts.nwk -p Mites_concat_AICc.nwk -o Lovette_Mites_concat_AICc.nex

python3 make_tangle.py -m amero_matrix.csv -s Lovette_hosts.nwk -p Mites_coi.nwk -o Lovette_Mites_coi.nex
```

I renamed the species in these and resaved these as: `Lovette_Mites_coi_jane_renamed.nex` and `Lovette_Mites_concatenated_jane_renamed.nex`
