#!/bin/bash

#SBATCH --job-name=00_PP_20210906_d_phylo
#SBATCH --partition=comp06
#SBATCH --output=00_PP_20210906_d_phylo_%j.txt
#SBATCH --error=00_PP_20210906_d_phylo_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00

module load java/sunjdk_1.8.0
module load fastqc/0.11.5
module load bbmap/38.82

module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7

cd /local_scratch/$SLURM_JOB_ID/

############ CONFIGURATION INFORMATION ############

############ LOCATIONS AND SAMPLE NAMES ############

## PROJECT DIRECTORY 
PROJECT_SAMPLE_DIR=/scrfs/storage/amatthews/20210816_projects

## ASSIGN 00_PP FOLDER TO QTRIM_DATA
QTRIM_DATA=${PROJECT_SAMPLE_DIR}/20210816_phylo/00_PP_20210902

## ADAPTER REMOVED
ADP_REM_BB=$QTRIM_DATA/Adapter_Removed_bb

## FASTQC DIR
A_Q_T_BB=$QTRIM_DATA/FASTQC/FASTQ_After_Q_Trim_bb



############ DEDUPE - one at a time... ############

############ MOWA_611_R_HC7TNCCX2_subsampled
dedupe.sh in1=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_R1.fastq in2=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_R2.fastq out=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_dd.fastq ac=false

# ac=false sets absorbContainment=false (which is necessary for paired reads)
# naming this output *_trim_dd so that my reformatted output can be *_trim_R1/R2 (to match endings with the other outputs that were not deduplicated)

reformat.sh in=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_bb_trim_dd.fastq out1=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_dd_bb_trim_R1.fastq out2=$ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_dd_bb_trim_R2.fastq

# so now when calling deduplicated samples, will need to add "_dd" to the end of their file names 

cd $A_Q_T_BB/MOWA_611_R_HC7TNCCX2

# run fastqc on these dd'd samples
fastqc $ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_dd_bb_trim_R1.fastq  $ADP_REM_BB/MOWA_611_R_HC7TNCCX2/MOWA_611_R_HC7TNCCX2_subsampled_dd_bb_trim_R2.fastq -o $PWD




############ BPWA_337_C_CAAGCCAA-CTTCGTTC
dedupe.sh in1=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_bb_trim_R1.fastq in2=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_bb_trim_R2.fastq out=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_bb_trim_dd.fastq ac=false

# ac=false sets absorbContainment=false (which is necessary for paired reads)
# naming this output *_trim_dd so that my reformatted output can be *_trim_R1/R2 (to match endings with the other outputs that were not deduplicated)

reformat.sh in=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_bb_trim_dd.fastq out1=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_dd_bb_trim_R1.fastq out2=$ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_dd_bb_trim_R2.fastq

# so now when calling deduplicated samples, will need to add "_dd" to the end of their file names 

cd $A_Q_T_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC

# run fastqc on these dd'd samples
fastqc $ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_dd_bb_trim_R1.fastq  $ADP_REM_BB/BPWA_337_C_CAAGCCAA-CTTCGTTC/BPWA_337_C_CAAGCCAA-CTTCGTTC_dd_bb_trim_R2.fastq -o $PWD




############ CAWA_366_A_CCGGAATA-CCTACTGA
dedupe.sh in1=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_bb_trim_R1.fastq in2=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_bb_trim_R2.fastq out=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_bb_trim_dd.fastq ac=false

# ac=false sets absorbContainment=false (which is necessary for paired reads)
# naming this output *_trim_dd so that my reformatted output can be *_trim_R1/R2 (to match endings with the other outputs that were not deduplicated)

reformat.sh in=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_bb_trim_dd.fastq out1=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_dd_bb_trim_R1.fastq out2=$ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_dd_bb_trim_R2.fastq

# so now when calling deduplicated samples, will need to add "_dd" to the end of their file names 

cd $A_Q_T_BB/CAWA_366_A_CCGGAATA-CCTACTGA

# run fastqc on these dd'd samples
fastqc $ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_dd_bb_trim_R1.fastq  $ADP_REM_BB/CAWA_366_A_CCGGAATA-CCTACTGA/CAWA_366_A_CCGGAATA-CCTACTGA_dd_bb_trim_R2.fastq -o $PWD




############ multiqc ############
# run multiqc on these dd'd samples (in addition to the others)
multiqc /scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/FASTQC/FASTQ_After_Q_Trim_bb





