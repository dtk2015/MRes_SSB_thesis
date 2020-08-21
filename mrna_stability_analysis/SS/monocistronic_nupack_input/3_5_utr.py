import pandas as pd
import numpy as np

#this code generates a csv of monocistornci UTR sequeces for input into nupakc

utr = pd.read_csv("UTR_5_3_sequence.csv")
utrcolnam = ['Operon Name',
 'mRNA_name',
 'Promoter Name',
 'TSS(+1)',
 'DNA_strand',
 'First_Gene',
 'Last_Gene',
 'Terminator Type [Terminator Position left, Terminator Position right]',
 'Coordinates',
 "Coordinates 5' UTR",
 "5' UTR Sequence",
 "Coordinates 3' UTR",
 "3' UTR Sequence"]

utr.columns = utrcolnam
#get the column of gene names for cds
#delete rows which have no information for 5' or 3' UTR (as want full length mRNA sequences)
utr.dropna(subset=["5' UTR Sequence"], inplace=True)
utr.dropna(subset=["3' UTR Sequence"], inplace=True)
#remove polycistronic genes (i.e. removed genes in dataset where first gene and last gene are not the same)
utrfin = utr.query("First_Gene == Last_Gene")
subopt = len(utrfin)*[1.5]
#add subopt value for getting suboptimal predicted secondary structures
utrfin["subopt"] = subopt
utrfin = utrfin.dropna(subset=['mRNA_name'])



utrfin["gene_name"] = utrfin["mRNA_name"]

ls = list(utrfin['mRNA_name'])

#this renames gene names to ensure genes with mutliple sequences are given separate names
for index,value in enumerate(ls):
    if value in ls[index+1:]:
        for new in range(0,200000000):
            if not  f'{new}{value}' in ls:
                ls[index] = f'{new}{value}'
                break

utrfin["mRNA_name"] = ls


#
UTR5 = utrfin[["mRNA_name", "5' UTR Sequence", 'gene_name',  "subopt"]]
UTR3 = utrfin[["mRNA_name", "3' UTR Sequence", 'gene_name', "subopt"]]


#move these csvs into the 3 and 5 folders and run csv_to_txt.py
UTR5.to_csv(r'5utr_seq.csv', index = False)
UTR3.to_csv(r'3utr_seq.csv', index = False)
