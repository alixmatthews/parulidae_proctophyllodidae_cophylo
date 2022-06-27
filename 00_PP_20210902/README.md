# Preprocessing of sequencing data for <i>Amerodectes</i> (et al.) phylogeny

Alix Matthews

Date: 20210902

## Computer and path of working directory
- AHPCC: ```/scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902```

---

### Initial Pre-processing 

- Check the quality of the data and remove poor quality data
-  `00_PP_20210902_a_phylo.slurm`
- AHPCC: ``` /scrfs/storage/amatthews/20210816_projects/20210816_phylo/slurms/00_PP_20210902_a_phylo.slurm ```

#### necessary files

- names of files that are to undergo initial pre-processing: `20210816_phylo_filenames.txt`


---

### BBDuk Pre-processing 

- Remove 10 bp from the 3' end (5 is for sure, but 5 more will probably be useful) - so cutting >140
- GC content is weird for CAWA (expected, only 2 mites in extraction), NAWA (expected, only 4 mites in extraction)
- Adapter content is highest for MOWA (expected, the Novogene sample)

- Does a second pass of pre-processing on the data
- `00_PP_20210902_b_phylo.slurm`
- AHPCC: `/scrfs/storage/amatthews/20210816_projects/20210816_phylo/slurms/00_PP_20210903_b_phylo.slurm`


#### necessary files

- names of files that are to undergo a second step of pre-processing: `20210816_phylo_filenames.txt`
- Adapters file necessary: `adapters.fa` 

---

### BBMap randomly subsample reads from MOWA Novogene sample 

- After 00_PP_20210903_b_phylo, I have 20.6 mil reads for MOWA (and between 4.9-8.4 M for other samples), so if I take 35% of these 20.6 M reads, I'll get around 7.21 M reads, which is in the range of my other samples
- I can use `reformat.sh` in bbmap to randomly subsample 35% of reads (while retaining pairs)

- `00_PP_20210903_c_snp.slurm`
- AHPCC: `/scrfs/storage/amatthews/20210816_projects/20210816_phylo/slurms/00_PP_20210903_c_snp.slurm`


---

### BBDuk deduplication 

- Remove duplicate sequences from certain samples
- `00_PP_20210906_d_snp.slurm`
- AHPCC: `/scrfs/storage/amatthews/20210816_projects/20210816_phylo/slurms/00_PP_20210906_d_snp.slurm`


---
### Need to manually change file structure for de-duplicated samples

- This is to make things easier when needing to reference only the deduplicated samples (in the names file). Unsure if I could have anticipated this and done so in one of the above scripts, but not going to dwell on it. Otherwise, the names files won't be able to reference/cd into the correct directory

- Manually did the following:
```
mkdir /scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/Adapter_Removed_bb/MOWA_611_R_HC7TNCCX2_subsampled_dd

# then move subsampled_dd data into this directory
# do not worry about moving non-dd'd subsampled data (it can stay in MOWA_611_R_HC7TNCCX2 directory)

mkdir /scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/Adapter_Removed_bb/BPWA_337_C_CAAGCCAA-CTTCGTTC_dd

# then move dd data into this directory

mkdir /scrfs/storage/amatthews/20210816_projects/20210816_phylo/00_PP_20210902/Adapter_Removed_bb/CAWA_366_A_CCGGAATA-CCTACTGA_dd

# then move dd data into this directory
```

---
### Output to proceed with
- names of files to proceed with, including those that have undergone deduplication (end in _dd): `20210816_phylo_filenames_post_dd.txt`
- These files are ready to move on: `$ADP_REM_BB/${value1}/${value1}_bb_trim_R1.fastq` and `$ADP_REM/${value1}/${value1}_bb_trim_R2.fastq`





