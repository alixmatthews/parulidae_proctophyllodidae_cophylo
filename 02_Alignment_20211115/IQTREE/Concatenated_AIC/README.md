## Files in this directory:
Various output files from IQTREE analysis using AIC as model selection criterion

- `Mites_tree_concat_AIC.best_scheme`: provides the best partitioning scheme found to avoid running ModelFinder again, coordinates for each partition of concatenated genes
- `Mites_tree_concat_AIC.best_scheme.nex`: provides the best partitioning scheme found to avoid running ModelFinder again, in .nexus format
- `Mites_tree_concat_AIC.bionj`: neighbor joining (BIONJ) tree in newick format
- `Mites_tree_concat_AIC.ckp.gz`: periodic IQ-TREE checkpoint file in case needed to resume an interrupted run
- `Mites_tree_concat_AIC.contree`: consensus tree constructed from all ultrafast bootstrap trees with ultrafast bootstrap support values 
- `Mites_tree_concat_AIC.iqtree`: IQ-TREE report
- `Mites_tree_concat_AIC.log`: printed screen log file of the entire run
- `Mites_tree_concat_AIC.mldist`: Likelihood distances from pairwise sequence comparisons based on model parameter estimates from an initial tree
- `Mites_tree_concat_AIC.model.gz`: log likelihoods for all models tested; checkpoint file to recover an interrupted model selection
- `Mites_tree_concat_AIC.splits.nex`: support values (%) for all bipartitions, computed as the occurence frequencies in the bootstrap trees
- `Mites_tree_concat_AIC.treefile`: best tree found by maximum-likelihood in newick format with ultrafast bootstrap support values 
