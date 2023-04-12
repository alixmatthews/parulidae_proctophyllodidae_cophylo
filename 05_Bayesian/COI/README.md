## Files in this directory:
Various input/output files for Bayesian analysis on COI dataset

- `04_PhyloBayes_chain1_coi_comp01.slurm`: first chain for PhyloBayes analysis; slurm file for HPC
- `04_PhyloBayes_chain1_coi_comp01_restart.slurm`: first chain, continuation, for PhyloBayes analysis; slurm file for HPC
- `04_PhyloBayes_chain2_coi_comp01.slurm`: second chain for PhyloBayes analysis; slurm file for HPC
- `04_PhyloBayes_chain2_coi_comp01_restart.slurm`: second chain, continuation, for PhyloBayes analysis; slurm file for HPC
- `chain12.bpdiff`: output from PhyloBayes after running bpcomp, displays values which are used to assess convergence of chains
- `chain12.bplist`: the list of weighted bipartitions from chains
- `chain12.con.tre`: the majority-rule posterior consensus tree for coi
- `coi_cat.fasta.mafftaligned_phylobayes.phy`: .fasta file of coi converted into a .phy file for input into PhyloBayes
