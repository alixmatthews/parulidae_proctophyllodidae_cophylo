'''
this will take names of genes from file and get the fasta file with that names. Concatenate fasta files coming from the same gene. 
'''

import os
import os.path
import glob
import sys



#file_dir = "/Users/aselawijeratne/Documents/Projects/scripts/python"

file_dir = sys.argv[1]
folder_pattern = sys.argv[2]
file_name = sys.argv[3]

#folder pattern for where the fasta files are:
#fileNames = glob.glob('**/*/*/*/*')

fileNames = glob.glob(folder_pattern)

#print fileNames

#sample_name = file_dir.split("/")[-1]

# file containing the names needs to be formatted.

#mt_genes = file_dir + "/" + "mt_genes.txt"

mt_genes = file_dir + "/" + file_name


with open(mt_genes) as input_files:
    for lines in input_files:
        cat_fasta = file_dir + "/" + lines.strip("\n") + "_cat.fasta"
        #print lines
        str_match = [s for s in fileNames if lines.strip("\n") in s]
        
        for seq in str_match:
            with open(seq) as input_fasta:
                for lines in input_fasta:
                    with open(cat_fasta, "a") as outFasta:
                        outFasta.write(lines)