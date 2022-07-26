#!/bin/bash

#SBATCH --job-name=04_PhyloBayes_chain4
#SBATCH --partition=comp72
#SBATCH --output=04_PhyloBayes_chain4_%j.txt
#SBATCH --error=04_PhyloBayes_chain4_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=20:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
conda activate phylobayes-mpi 


# I did this (below) to create a PhyloBayes environment in the home login node beforehand
# module load python/anaconda-3.8
# source /share/apps/bin/conda-3.8.sh
# mamba create -n phylobayes-mpi phylobayes-mpi 

# PhyloBayes is version 1.9


cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"/DIRECTORIES

## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects


### ASSIGN PHYLO_DIR
PHYLO_DIR=${PROJECT_DIR}/20210816_phylo/03_PhyloBayes


### DATAFILE 
DATAFILE=${PHYLO_DIR}/Mites_concat_phylobayes.phy


### CHAIN NAME
CHAIN_NAME=chain4




##----------------------------------------------------------------------------------------
## Run PhyloBayes

#### -np: number of parallel processes (32 for AHPCC)
#### -d: datafile
#### -cat -dp: infinite mixture (Dirichlet process) of equilibrium freq profiles
#### -gtr: GTR model
#### -dagm: discrete gamma (4 is default)
#### -s: save all

cd ${PHYLO_DIR}

mpirun -np 32 pb_mpi -s -d ${DATAFILE} -cat -dp -gtr -dgam 4 ${CHAIN_NAME}




