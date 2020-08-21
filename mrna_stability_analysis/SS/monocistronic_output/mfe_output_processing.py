#this code creates a dataframe from  the mfe output files from nupack in '3' and '5' folder
#NOTE: these mfe files could not be uploaded to GitHub as took up too much space
#therefore this code will not run wihtout first creating these mfe files
#the output of this code are '3_UTR_nupack.csv' and '5_UTR_nupack.csv'

import pandas as pd 
import os 

#this code processes nupack output files into a csv

##############################################################################
##################################5utr########################################
##############################################################################

#########
###mfe###
#########

#get list of gene names from titles of files
five_mfe_gene_name = os.listdir(r"C:\Users\danyn\OneDrive\Desktop\Github\mrna_stability_analysis\SS\monocistronic_output\5\mfe") #path to location of 5' utr nupakc output files

#specify file paths
five_mfe_folderpath = r"C:\Users\danyn\OneDrive\Desktop\Github\mrna_stability_analysis\SS\monocistronic_output\5\mfe" #path to location of 5' nupakc output files
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

five_mfe_df3.to_csv("5UTR_nupack.csv", mode = 'w', index=False)

##############################################################################
##################################3utr########################################
##############################################################################
#get list of gene names from titles of files
three_mfe_gene_name = os.listdir(r"C:\Users\danyn\OneDrive\Desktop\Github\mrna_stability_analysis\SS\monocistronic_output\3\mfe") #path to location 3' utr of nupakc output files

#specify file paths
three_mfe_folderpath = r"C:\Users\danyn\OneDrive\Desktop\Github\mrna_stability_analysis\SS\monocistronic_output\3\mfe" #path to location of 3' utr nupakc output files
three_mfe_filepaths  = [os.path.join(three_mfe_folderpath, name) for name in os.listdir(three_mfe_folderpath)]
three_mfe_all_files = []

#create a lsit of each line in each file
for path in three_mfe_filepaths:
    with open(path, 'r') as f:
        file = f.readlines()
        three_mfe_all_files.append(file)

#create a df of the list of lists
three_mfe_df=pd.DataFrame(three_mfe_all_files)

#create a df of the gene names
three_mfe_df1=pd.DataFrame(three_mfe_gene_name)

#create a df from rows containign sequence, length, mfe, and structure
three_mfe_df2 = three_mfe_df[[4,14,15,16]].copy()

#create a df from df2 and gene name df
three_mfe_df3 = pd.concat([three_mfe_df1, three_mfe_df2], axis = 1)

#remove unwanted strings from columns
three_mfe_df3 = three_mfe_df3.replace({"%":""}, regex=True)
three_mfe_df3 = three_mfe_df3.replace({" ":""}, regex=True)
three_mfe_df3 = three_mfe_df3.replace(r'\n','', regex=True) 
three_mfe_df3[0] = three_mfe_df3[0].str.replace(r'3','')
three_mfe_df3[0] = three_mfe_df3[0].str.replace(r'.mfe','')
three_mfe_df3[4] = three_mfe_df3[4].str.replace(r"Sequence:",'')

#convert string numbers to integers
three_mfe_df3 = three_mfe_df3.apply(pd.to_numeric, errors = "ignore")

#add column name
three_mfe_df3.columns = ["gene_name", "sequence", "length", "mfe", "mfe_secondary_structure"]

three_mfe_df3.to_csv("3UTR_nupack.csv", mode = 'w', index=False)