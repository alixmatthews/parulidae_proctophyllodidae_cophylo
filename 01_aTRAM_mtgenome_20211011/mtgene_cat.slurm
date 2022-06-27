#!/bin/bash

#SBATCH --job-name=mtgene_cat
#SBATCH --partition=comp01
#SBATCH --output=mtgene_cat_%j.txt
#SBATCH --error=mtgene_cat_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:09:00

module load python/2.7.15

cd /local_scratch/$SLURM_JOB_ID/

############################################################
############ CONFIGURATION INFORMATION ############
############ DIRECTORIES ############
############################################################


## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_phylo/01_aTRAM_mtgenome_20211011

cd ${PROJECT_DIR}

python gene_cat.py "/scrfs/storage/amatthews/20210816_projects/20210816_phylo/01_aTRAM_mtgenome_20211011" "/scrfs/storage/amatthews/20210816_projects/20210816_phylo/01_aTRAM_mtgenome_20211011/*_out/*_mtgenome/stitch/*" "mt_genes.txt"





