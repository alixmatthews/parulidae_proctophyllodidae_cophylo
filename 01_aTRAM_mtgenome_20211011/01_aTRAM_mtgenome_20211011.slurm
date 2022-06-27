#!/bin/bash

#SBATCH --job-name=01_aTRAM_mtgenome_20211011
#SBATCH --partition=comp72
#SBATCH --output=01_aTRAM_mtgenome_20211011_%j.txt
#SBATCH --error=01_aTRAM_mtgenome_20211011_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=36:00:00

module load exonerate/2.4.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate aTRAM

cd /local_scratch/$SLURM_JOB_ID/

############################################################
############ CONFIGURATION INFORMATION ############
############ LOCATIONS AND SAMPLE NAMES ############
############################################################

## SAMPLE FILE NAMES 
# this file has all the Kays samples that have been pre-processed up until this date (i.e., the CAWA and BPWA still may need some PP post deduplicating, but just giving it a shot now...)
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_phylo/20210816_phylo_filenames_post_dd.txt" 

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_phylo

## ATRAM DIR
## mkdir beforehand
## contains mtgenome reference (Proctophyllodes_miliariae)
ATRAM_DIR=$PROJECT_DIR/01_aTRAM_mtgenome_20211011

## ATRAM DATABASES
ATRAM_DB=$ATRAM_DIR/aTRAM_db
mkdir ATRAM_DB


############ SET UP PER SAMPLE DIRS ############

while IFS=" " read -r value1
do {


## FASTQ FILE DIRECTORY
SAMPLE_FASTQ=$PROJECT_DIR/00_PP_20210902/Adapter_Removed_bb/${value1}

## ATRAM_DB DIR
## this line actually doesn't do anything meaningful (forgot to put a $ before {}), but keeping because the rest of the script works and don't want removing this to mess it up!
mkdir $ATRAM_DB/{value1}

## SAMPLE OUTPUT / mtgenome / (ALL, FILTERED, STITCH)
SAMPLE_OUT=$ATRAM_DIR/${value1}_out
mkdir $SAMPLE_OUT
SAMPLE_OUT_mtgenome=$SAMPLE_OUT/${value1}_mtgenome
mkdir -p $SAMPLE_OUT_mtgenome
mkdir $SAMPLE_OUT_mtgenome/all
mkdir $SAMPLE_OUT_mtgenome/filtered
mkdir $SAMPLE_OUT_mtgenome/stitch



############################################################
############ RUN aTRAM ############
############################################################

############ PREPROCESSOR ############

atram_preprocessor.py --blast-db=$ATRAM_DB/${value1}/${value1} --end-1=$SAMPLE_FASTQ/${value1}_bb_trim_R1.fastq  --end-2=$SAMPLE_FASTQ/${value1}_bb_trim_R2.fastq
## second /${value1} on --blast-db line is referring to the prefix to use for the file names it creates in this step


############ ASSEMBLER ############
############ mtgenome ############

atram.py --blast-db=$ATRAM_DB/${value1}/${value1} --query-split=$ATRAM_DIR/mt_ref.fasta --assembler=spades --protein --db-gencode=5 --output-prefix=$SAMPLE_OUT_mtgenome/all/
## second /${value1} on --blast-db line is referring to the prefix to use/search for the file names it created in the previous step

cd $SAMPLE_OUT_mtgenome
find . -name "*filtered_contigs.fasta" -exec mv -i {} -t "./filtered/" ";"
## move filtered_contigs files to a new folder called "filtered" (separate them from non-filtered)



############ STITCHER ############
############ mtgenome ############

atram_stitcher.py --assemblies-dir=$SAMPLE_OUT_mtgenome/filtered --reference-genes=$ATRAM_DIR/mt_ref.fasta --taxa=/scrfs/storage/amatthews/20210816_projects/20210816_phylo/20210816_phylo_filenames_post_dd.txt --output-prefix=$SAMPLE_OUT_mtgenome/stitch/${value1}



############ FINISH aTRAM ############

} done <"$filenames"

