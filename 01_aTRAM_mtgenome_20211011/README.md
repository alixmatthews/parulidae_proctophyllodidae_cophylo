# Order of operations post 00-PP

## Gene assembly from WGS
`01_aTRAM_mtgenome_20211011.slurm` - Assemble full mtgenome (13 PCGs) for all 33 Kays samples (*Amerodectes* and *Tyrannidectes*), using mitochondrial reference is amino acid translated mt genes from *Proctophyllodes miliariae*, and aTRAM

`mt_genes.txt` - file containing the names of the 13 mito PGCs

`gene_cat.py` - python script to extract the fasta files by gene from each samples subdirectories and then concatenate them together (headers in the output geneID_cat.fasta file are the sample names). This file is located directly in the $PROJECT_DIR as designated in `mtgene_cat.slurm` below.

`mtgene_cat.slurm` - slurm script to submit `gene_cat.py` as a job on AHPCC.  

The end result is several .fasta files named `$gene_cat.fasta` and these are ready to move on to 02_Alignment!
