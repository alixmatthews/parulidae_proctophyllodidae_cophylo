#!/bin/bash

#SBATCH --job-name=00_PP_20210902_a_phylo
#SBATCH --partition=comp72
#SBATCH --output=00_PP_20210902_a_phylo_%j.txt
#SBATCH --error=00_PP_20210902_a_phylo_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00

module load java/sunjdk_1.8.0
module load fastqc/0.11.5

module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7

cd /local_scratch/$SLURM_JOB_ID/

############ CONFIGURATION INFORMATION ############

############ LOCATIONS AND SAMPLE NAMES ############

## SAMPLE FILE NAMES 
## value1 is sample name
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_phylo/20210816_phylo_filenames.txt"

## PROJECT DIRECTORY 
PROJECT_SAMPLE_DIR=/scrfs/storage/amatthews/20210816_projects

## WHERE THE RAW DATA ARE
RAW_DATA=${PROJECT_SAMPLE_DIR}/raw_data

## ASSIGN 00_PP FOLDER TO QTRIM_DATA
QTRIM_DATA=${PROJECT_SAMPLE_DIR}/20210816_phylo/00_PP_20210902

## MAKE NEW DIRECTORIES
mkdir -p $QTRIM_DATA ## necessary to mk multiple dirs within QTRIM_DATA
mkdir $QTRIM_DATA/FASTQC ## FASTQC
mkdir $QTRIM_DATA/FASTQC/Before_Q_Trim ## mkdir within FASTQC called Before_Q_Trim

## ADAPTER REMOVED
ADP_REM=$QTRIM_DATA/Adapter_Removed 
mkdir $ADP_REM ## mkdir to the ADP_REM object 

## FASTQC AFTER QUALITY TRIM
A_Q_T=$QTRIM_DATA/FASTQC/FASTQ_After_Q_Trim ## make an object called A_Q_T within the $QTRIM_DATA/FASTQC directory which will be a folder named FASTQ_After_Q_Trim
mkdir $A_Q_T ## mkdir to the A_Q_T object


############ SOFTWARE SETTINGS ############

## QUALITY STUFFS
FASTQ_TYPE="sanger"
QUALIY_THRESHOLD=30
MIN_LENGTH=75 ## can optimize this later if desired
MIN_OVERLAP=3 ## default is 3, but may need to change if 6 if lots of seqs get thrown out
# TRIM_R1=-5 ## trim last (negative) 5 bases unambiguously from R1
# TRIM_R2=-5 ## trim last (negative) 5 bases unambiguously from R2

## ADAPTER INFO THAT NEEDS TO BE REMOVED
ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG  ## i7 TruSeq adapter (for forward read) (NNNNNNNN= 8 nt index)
U_ADAPTER_SEQ=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNAGATCTCGGTGGTCGCCGTATCAT ## i5 TruSeq universal Illumina adapter (for reverse read)

############ START THE LOOP THROUGH FILES ############

while IFS=" " read -r value1
do {


############ ASSIGN FILES TO OBJECTS ############

## each line of the "file" is read with the ${line} syntax

FIRST_SAMPLE_LOC=${RAW_DATA}/${value1}_L002_R1_001.fastq.gz
SECOND_SAMPLE_LOC=${RAW_DATA}/${value1}_L002_R2_001.fastq.gz


############ FASTQC ############

mkdir $QTRIM_DATA/FASTQC/Before_Q_Trim/${value1} ## mkdir of a folder called whatever my sample name is (${value1}) within the Before_Q_Trim folder
cd $QTRIM_DATA/FASTQC/Before_Q_Trim/${value1} ## cd into that sample's folder

## run fastqc and put the output in the pwd ($QTRIM_DATA/FASTQC/Before_Q_Trim/${line})
fastqc $FIRST_SAMPLE_LOC $SECOND_SAMPLE_LOC -o $PWD


############ REMOVE ADAPTERS ############

mkdir $ADP_REM/${value1} ## mkdir of a folder called whatever my sample name is (${value1}) within the $QTRIM_DATA/Adapter_Removed folder
cd $ADP_REM/${value1} ## cd into that sample's folder

## run cutadapt and put the output in the pwd ($ADP_REM/${value1})
cutadapt -q $QUALIY_THRESHOLD -m $MIN_LENGTH -a $ADAPTER_SEQ -A $U_ADAPTER_SEQ -O $MIN_OVERLAP -o ${value1}_trim_R1.fastq -p ${value1}_trim_R2.fastq $FIRST_SAMPLE_LOC $SECOND_SAMPLE_LOC 
# can add -u $TRIM_R1 -U $TRIM_R2 if choose to adjust these software settings above


############ FASTQC AFTER TRIMMING ############

mkdir $A_Q_T/${value1} ## mkdir of a folder called whatever my sample name is (${value1}) within the $QTRIM_DATA/FASTQC/FASTQ_After_Q_Trim folder
cd $A_Q_T/${value1} ## cd into that sample's folder

fastqc $ADP_REM/${value1}/${value1}_trim_R1.fastq $ADP_REM/${value1}/${value1}_trim_R2.fastq -o $PWD

} done <"$filenames"

############ MULTIQC ############

## run multiqc - just give it the path to fastqc files
multiqc /scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/FASTQC/FASTQ_After_Q_Trim
