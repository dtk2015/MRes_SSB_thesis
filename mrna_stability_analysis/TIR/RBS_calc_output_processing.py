import pandas as pd
import re

#This code generates whole dataframe of all TIR from RBS calculator subsets

######################################################################
######################################################################
######################################################################
#concatenate rbs sub df into single df
#WARNING ONLY RUN ONCE AS TAKES ALOT OF MEMORY AND TIME TO DO THIS
#open the 5 rbs clacualtor csvs
rbs0 = pd.read_csv("RBS_Calculator_Predictions_sub0 result total.csv")
rbs1 = pd.read_csv("RBS_Calculator_Predictions_sub1 result total.csv")
rbs2 = pd.read_csv("RBS_Calculator_Predictions_sub2 result total.csv")
rbs3 = pd.read_csv("RBS_Calculator_Predictions_sub3 result total.csv")
rbs4 = pd.read_csv("RBS_Calculator_Predictions_sub4 result total.csv")
#concatenate them togetehr
rbs = pd.concat([rbs0, rbs1, rbs2, rbs3, rbs4])


#filter values
filter_vals=["TIR","ORF", "dG_total", "dG_mRNA_rRNA", "dG_mRNA", "dG_spacing"]

#rename mRNA_name column to name
rbs=rbs.rename(columns = {'mRNA_name':'name'})

dfs = {
    filter_name: rbs.set_index(["name", "mRNA_sequence"]).filter(like=filter_name)
    for filter_name in filter_vals}

tir = dfs['TIR']
orf = dfs['ORF']
dG_mRNA_rRNA = dfs['dG_mRNA_rRNA']
dG_mRNA = dfs['dG_mRNA']
dG_spacing = dfs['dG_spacing']



# #remove non-digits from postion columns
# #create list of columns in df
col_name_ls = list(tir)
# #create new lsit to appedn digit only column names (this list laready has first two column names)
col_name_ls_2 = []

#remove non-digits from these column names
for i in col_name_ls:
    no_dig = re.sub('[^0-9]','', i)
    col_name_ls_2.append(no_dig)
 #change column names to no digit column names
tir.columns = col_name_ls_2

#open original df we got utr sequecnes from
ecoli_df =pd.read_csv("Supplementary_Data_1_Ecoli_features_and_measurements.csv") 
ecoli_df_2 = ecoli_df
#remove spaces
ecoli_df['utr_cds_seq'] = ecoli_df['utr_cds_seq'].str.strip()
#check position start codon starts at
ecoli_df['start_codon'] = [utr_cds_seq.index(cds) for utr_cds_seq, cds in zip(ecoli_df['utr_cds_seq'], ecoli_df['cds'])]
#check to verify this is a start codon
ecoli_df["start_codon_str"] = ecoli_df['utr_cds_seq'].apply(lambda x: x[51:54])
print(list(ecoli_df["start_codon_str"]))

#merge TIR df with ecoli info df
full_df = tir.merge(ecoli_df,on='name',how='left') 

# #look up start position of cds in UTR sequence
#full_df['start_codon'] = [mRNA_sequence.index(cds) for mRNA_sequence, cds in zip(rbs['mRNA_sequence'], full_df['cds'])]
#print(full_df["start_codon"])
full_df['TIR']=full_df.lookup(full_df.index,full_df['start_codon'].astype(str))

#create new df of jsut the TIR and name
tir_df = full_df[["name", "TIR"]].copy()
#merge TIR df with ecoli info df
full_df = ecoli_df_2.merge(tir_df,on='name',how='left') 

full_df.to_csv('RBS_calc.csv')