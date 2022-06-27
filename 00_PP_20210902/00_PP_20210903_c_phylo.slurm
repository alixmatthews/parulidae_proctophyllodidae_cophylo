#!/bin/bash

#SBATCH --job-name=00_PP_20210903_c_phylo
#SBATCH --partition=comp01
#SBATCH --output=00_PP_20210903_c_phylo_%j.txt
#SBATCH --error=00_PP_20210903_c_phylo_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00

module load java/sunjdk_1.8.0
module load fastqc/0.11.5
module load bbmap/38.82

module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7

cd /local_scratch/$SLURM_JOB_ID/

## SAMPLE DIRECTORY 
PROJECT_SAMPLE_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/Adapter_Removed_bb/MOWA_611_R_HC7TNCCX2


reformat.sh in1=${PROJECT_SAMPLE_DIR}/MOWA_611_R_HC7TNCCX2_bb_trim_R1.fastq in2=${PROJECT_SAMPLE_DIR}/MOWA_611_R_HC7TNCCX2_bb_trim_R2.fastq out1=${PROJECT_SAMPLE_DIR}/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_R1.fastq out2=${PROJECT_SAMPLE_DIR}/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_R2.fastq samplerate=0.35

