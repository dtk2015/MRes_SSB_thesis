#this code creates a dataframe from  the mfe output files from nupack in 'mfe' folder
#NOTE: these mfe files could not be uploaded to GitHub as took up too much space
#therefore this code will not run wihtout first creating these mfe files
#the output of this code is 'kim_2008_nupack_processed.csv'
import pandas as pd 
import os 
##############################################################################
##################################5utr########################################
##############################################################################

#########
###mfe###
#########

#get list of gene names from titles of files
five_mfe_gene_name = os.listdir(r"D:\Imperial\MRes_project2\SS\kim2008_output\mfe") #directory of neupack outptu files

#specify file paths
five_mfe_folderpath = r"D:\Imperial\MRes_project2\SS\kim2008_output\mfe"  #directory of neupack outptu files
five_mfe_filepaths  = [os.path.join(five_mfe_folderpath, name) for name in os.listdir(five_mfe_folderpath)]
five_mfe_all_files = []

#create a lsit of each line in each file
for path in five_mfe_filepaths:
    with open(path, 'r') as f:
        file = f.readlines()
        five_mfe_all_files.append(file)

#create a df of the list of lists
five_mfe_df=pd.DataFrame(five_mfe_all_files)

#create a df of the gene names
five_mfe_df1=pd.DataFrame(five_mfe_gene_name)

#create a df from rows containign sequence, length, mfe, and structure
five_mfe_df2 = five_mfe_df[[4,13,14,15]].copy()

#create a df from df2 and gene name df
five_mfe_df3 = pd.concat([five_mfe_df1, five_mfe_df2], axis = 1)

#remove unwanted strings from columns
five_mfe_df3 = five_mfe_df3.replace({"%":""}, regex=True)
five_mfe_df3 = five_mfe_df3.replace({" ":""}, regex=True)
five_mfe_df3 = five_mfe_df3.replace(r'\n','', regex=True) 
five_mfe_df3[0] = five_mfe_df3[0].str.replace(r'5','')
five_mfe_df3[0] = five_mfe_df3[0].str.replace(r'.mfe','')
five_mfe_df3[4] = five_mfe_df3[4].str.replace(r"Sequence:",'')

#convert string numbers to integers
five_mfe_df3 = five_mfe_df3.apply(pd.to_numeric, errors = "ignore")

#add column name
five_mfe_df3.columns = ["gene_name", "sequence", "length", "mfe", "mfe_secondary_structure"]


df = pd.read_csv('kim_2008_seq.csv')
df.columns = ['TSS position',
 'Strand',
 'Gene ID',
 'Gene name',
 'Gene start',
 'Gene end',
 'sequence',
 'seq_len',
 'gene name',
 'subopt']

x1 = five_mfe_df3.merge(df,on='sequence',how='left') 

x2 = x1[['gene name', 'sequence', 'length', 'mfe', 'mfe_secondary_structure']]

#add column name
x2.columns = ["gene_name", "sequence", "length", "mfe", "mfe_secondary_structure"]
x2.sort_values(by=['gene_name'], inplace=True)

#save output csv
x2.to_csv("5UTR_nupack.csv", mode = 'w', index=False)
