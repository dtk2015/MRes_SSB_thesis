import pandas as pd
import numpy as np

#PSA
circ_PSA_scaled_df = pd.read_csv (r'PSA_scaled.txt', sep='\t') 
circ_PSA_scaled_df.columns = ['species', 'k1', 'kd1', 'kt3', 'kt5', 'k2', 'kd2']
circ_PSA_scaled_df['species'].replace('[l1]', 'l1RNA',inplace=True)
circ_PSA_scaled_df['species'].replace('[l13]', 'l13RNA',inplace=True)
circ_PSA_scaled_df['species'].replace('[l15]', 'l15RNA',inplace=True)
circ_PSA_scaled_df['species'].replace('[l135]', 'l135RNA',inplace=True)
circ_PSA_scaled_df['species'].replace('[c]', 'cmRNA',inplace=True)
circ_PSA_scaled_df =  circ_PSA_scaled_df.set_index('species')

lin_PSA_scaled_df = pd.read_csv (r'lin_PSA_scaled.txt', sep='\t') 
lin_PSA_scaled_df.columns = ['species', 'k1', 'skd1', 'kt3', 'kt5', 'k2', 'kd2', 'sk1', 'kd1', 'kc', 'kd', 'k3']
lin_PSA_scaled_df['species'].replace('[l1]', 'l1RNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[l13]', 'l13RNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[l15]', 'l15RNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[l135]', 'l135RNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[c]', 'cmRNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[s]', 'siRNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[cs]', 'cmRNA:siRNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[l2]', 'l2RNA',inplace=True)
lin_PSA_scaled_df['species'].replace('[sa]', 'siRNAa',inplace=True)
lin_PSA_scaled_df['species'].replace('[sb]', 'siRNAb',inplace=True)
lin_PSA_scaled_df =  lin_PSA_scaled_df.set_index('species')
lin_PSA_scaled_df_1 = lin_PSA_scaled_df[['sk1', 'kd1','k3']]
lin_PSA_scaled_df_1 = lin_PSA_scaled_df_1.drop(['l1RNA', 'l13RNA', 'l15RNA', 'l135RNA', 'siRNAa', 'siRNAb'], axis=0)
#create heatmap 
import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10,10))
sns.heatmap(circ_PSA_scaled_df,ax=ax, cmap="rocket", annot=True, fmt='.2g')
ax.set_ylabel('')    
ax.set_xlabel('')
plt.show()


fig, ax = plt.subplots(figsize=(10,10))
ax = sns.heatmap(lin_PSA_scaled_df_1,ax=ax, cmap="rocket", annot=True)
ax.set_ylabel('')    
ax.set_xlabel('')
plt.show()


















