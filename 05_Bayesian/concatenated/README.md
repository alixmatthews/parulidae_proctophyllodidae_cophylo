## Files in this directory:
Various input/output files for Bayesian analysis on concatenated dataset

- `04_PhyloBayes_chain1.slurm`: first chain for PhyloBayes analysis; slurm file for HPC
- `04_PhyloBayes_chain4.slurm`: second chain for PhyloBayes analysis; slurm file for HPC
- `Mites_concat_phylobayes.phy`: .fasta file of concatenated genes converted into a .phy file for input into PhyloBayes
- `chain14.bpdiff`: output from PhyloBayes after running bpcomp, displays values which are used to assess convergence of chains
- `chain14.bplist`: the list of weighted bipartitions from chains
- `chain14.con.tre`: the majority-rule posterior consensus tree for concatenated dataset
