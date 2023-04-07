# Using aTRAM to assemble mt genes after pre-processing steps

## Files in this directory:

- `01_aTRAM_mtgenome_20211011.slurm`
- `gene_cat.py`
- `mt_genes.txt`
- `mt_ref.fasta`
- `mtgene_cat.slurm`


---

### Info about files


`01_aTRAM_mtgenome_20211011.slurm` - Assemble full mtgenome (13 PCGs) for all 33 Kays samples (*Amerodectes* and *Tyrannidectes*) using aTRAM... the mitochondrial reference (`mt_ref.fasta`) is the amino acid translated mt genes from *Proctophyllodes miliariae*

`mt_genes.txt` - file containing the names of the 13 mito PGCs

`gene_cat.py` - python script to extract the fasta files by gene from each samples subdirectories and then concatenate them together (headers in the output geneID_cat.fasta file are the sample names). This file is located directly in the $PROJECT_DIR as designated in `mtgene_cat.slurm` below.

`mtgene_cat.slurm` - slurm script to submit `gene_cat.py` as a job on AHPCC.  

The end result is several .fasta files named `$gene_cat.fasta` and these are ready to move on to 02_Alignment! Can be found in `RESULTS` directory here.
