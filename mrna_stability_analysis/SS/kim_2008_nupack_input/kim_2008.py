import pandas as pd

#this code generates a csv of 5' UTR sequeces from the kim 2008 study for input into nupakc


#open genome file as string NOTE I edited the txt genome file to remove the identifier manually 
#WARNING DO NOT PRINT SEQUENCE AS MY COMPUTER IS SHIT AND IT TAKES AGES TO PRINT OUT FULL GENOME
with open('Genome.txt', 'r') as file:
    genome = file.read().replace('\n', '')

#fucntion to get reverse compelemtn
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join([complement[base] for base in dna[::-1]])


#open dataset with UTR cooridnates
utr = pd.read_csv("kim_2008.csv")

#deete unecessary columns
del utr['Chromosome']
del utr['TSS reads']
del utr['Gene product']
del utr['COG function']
del utr["5' UTR length"]
del utr['Sequence -50 nt upstream + TSS (51 nt)']

#create empty dataframe for reverse strand mRNAs
utrrev = pd.DataFrame(columns=utr.columns)
#populate empty df with - strand mRNA info
cond = utr.Strand == "-"
rows = utr.loc[cond, :]
utrrev = utrrev.append(rows, ignore_index=True)
utr.drop(rows.index, inplace=True)

#get sequence of UTR from coordinates and add 25nt into cds, get length of sequence
utr['seq_UTR+25cds'] = [genome[a-1:b-1+25] for a,b in zip(utr["TSS position"], utr["Gene start"])]
utr["seq_len"] = utr["seq_UTR+25cds"].str.len()
#get sequence of UTR from coordinates and add 25nt into cds, get length of sequence and reverse complement seq (as - strand)
utrrev['seq_UTR+25cds'] = [genome[a-1-25:b-1] for a,b in zip(utrrev["Gene end"], utrrev["TSS position"])]
utrrev["seq_len"] = utrrev["seq_UTR+25cds"].str.len()
utrrev['seq_UTR+25cds'] = utrrev['seq_UTR+25cds'].apply(reverse_complement)

#concatenate two df together
utr = pd.concat([utr, utrrev])

utr["gene name"] = utr["Gene ID"]

ls = list(utr['Gene ID'])

#ensures genes with mutliple UTR sequences have unqiue names
for index,value in enumerate(ls):
    if value in ls[index+1:]:
        for new in range(0,200000000):
            if not  f'{new}{value}' in ls:
                ls[index] = f'{new}{value}'
                break

#change gene ids
utr["Gene ID"] = ls

#add subopt column
utr['subopt'] = 1.5
#MOVE THIS TO 'kim_2008' and run csv_to_tx.py
utr.to_csv(r'kim_2008_seq.csv', index = False)