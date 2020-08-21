import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

#sampling
samp_circ_df = pd.read_csv (r'whole_sample_circ_10.txt', sep='\t')
del samp_circ_df['Time']
samp_lin_df = pd.read_csv (r'whole_sample_lin_10.txt', sep='\t')
del samp_lin_df['Time']

#sum to egt total RNA in system
samp_circ_df['Total RNA'] = samp_circ_df[[ '[l1]','[l13]','[l15]','[l135]','[c]']].sum(axis=1)
samp_lin_df['Total RNA'] = samp_lin_df[[ '[l1]','[l13]','[l15]','[l135]','[c]','[cs]','[l2]',]].sum(axis=1)

#cirularisation
#get range of each rate constant
samp_circ_max_k1 = samp_circ_df["Values[k1].InitialValue"].max()
samp_circ_min_k1 = samp_circ_df["Values[k1].InitialValue"].min()
k1_range = samp_circ_max_k1 - samp_circ_min_k1

samp_circ_max_kd1 = samp_circ_df["Values[kd1].InitialValue"].max()
samp_circ_min_kd1 = samp_circ_df["Values[kd1].InitialValue"].min()
kd1_range = samp_circ_max_kd1 - samp_circ_min_kd1

samp_circ_max_kt3 = samp_circ_df["Values[kt3].InitialValue"].max()
samp_circ_min_kt3 = samp_circ_df["Values[kt3].InitialValue"].min()
kt3_range = samp_circ_max_kt3 - samp_circ_min_kt3

samp_circ_max_kt5 = samp_circ_df["Values[kt5].InitialValue"].max()
samp_circ_min_kt5 = samp_circ_df["Values[kt5].InitialValue"].min()
kt5_range = samp_circ_max_kt5 - samp_circ_min_kt5

samp_circ_max_k2 = samp_circ_df["Values[k2].InitialValue"].max()
samp_circ_min_k2 = samp_circ_df["Values[k2].InitialValue"].min()
k2_range = samp_circ_max_k2 - samp_circ_min_k2

samp_circ_max_kd2 = samp_circ_df["Values[kd2].InitialValue"].max()
samp_circ_min_kd2 = samp_circ_df["Values[kd2].InitialValue"].min()
kd2_range = samp_circ_max_kd2 - samp_circ_min_kd2

samp_circ_df = samp_circ_df.sort_values(by=['[c]'])

#create df for normalsied rates
column_names = ["normalised k1"]
samp_circ_rates_df = pd.DataFrame(columns = column_names)

#subtract the lower bound of the range from eahc rate in the df (normalisaiton)
samp_circ_rates_df["normalised k1"] = samp_circ_df["Values[k1].InitialValue"]-samp_circ_min_k1
samp_circ_rates_df["normalised kd1"] = samp_circ_df["Values[kd1].InitialValue"]-samp_circ_min_kd1
samp_circ_rates_df["normalised kt3"] = samp_circ_df["Values[kt3].InitialValue"]-samp_circ_min_kt5
samp_circ_rates_df["normalised kt5"] = samp_circ_df["Values[kt5].InitialValue"]-samp_circ_min_kt3
samp_circ_rates_df["normalised k2"] = samp_circ_df["Values[k2].InitialValue"]-samp_circ_min_k2
samp_circ_rates_df["normalised kd2"] = samp_circ_df["Values[kd2].InitialValue"]-samp_circ_min_kd2

#divide each rate cosntant column by range of rate constant (normalisation)
samp_circ_rates_df["normalised k1"] = samp_circ_rates_df["normalised k1"].div(k1_range)
samp_circ_rates_df["normalised kd1"] = samp_circ_rates_df["normalised kd1"].div(kd1_range)
samp_circ_rates_df["normalised kt3"] = samp_circ_rates_df["normalised kt3"].div(kt3_range)
samp_circ_rates_df["normalised kt5"] = samp_circ_rates_df["normalised kt5"].div(kt5_range)
samp_circ_rates_df["normalised k2"] = samp_circ_rates_df["normalised k2"].div(k2_range)
samp_circ_rates_df["normalised kd2"] = samp_circ_rates_df["normalised kd2"].div(kd2_range)

#creeate df of just species steady states
samp_circ_species_df = samp_circ_df[['[l1]','[l13]','[l15]','[l135]','[c]', 'Total RNA']]
#rename columns
samp_circ_species_df.columns = ['l1RNA', 'l13RNA', 'l15RNA', 'l135RNA', 'cmRNA', 'Total RNA']

#create heatmap of normalised rates and species steady states
import seaborn as sns
fig, (ax,ax2) = plt.subplots(ncols=2, figsize=(5,10))
fig.subplots_adjust(wspace=0.02)
sns.heatmap(samp_circ_rates_df, cmap="rocket", ax=ax, cbar=False, yticklabels=False)
fig.colorbar(ax.collections[0], ax=ax,location="left", use_gridspec=False, pad=0.2)
sns.heatmap(samp_circ_species_df, cmap="Blues", ax=ax2, cbar=False, yticklabels=False)
fig.colorbar(ax2.collections[0], ax=ax2,location="right", use_gridspec=False, pad=0.2)
ax2.yaxis.tick_right()
ax.tick_params(left=False, bottom=False)
ax2.tick_params(left=False, bottom=False)
plt.show()

#linearisation
#get range of each rate constant
samp_lin_max_lk1 = samp_lin_df["Values[lk1].InitialValue"].max()
samp_lin_min_lk1 = samp_lin_df["Values[lk1].InitialValue"].min()
lk1_lin_range = samp_lin_max_lk1 - samp_lin_min_lk1

samp_lin_max_kd1 = samp_lin_df["Values[kd1].InitialValue"].max()
samp_lin_min_kd1 = samp_lin_df["Values[kd1].InitialValue"].min()
kd1_lin_range = samp_lin_max_kd1 - samp_lin_min_kd1

samp_lin_max_k3 = samp_lin_df["Values[k3].InitialValue"].max()
samp_lin_min_k3 = samp_lin_df["Values[k3].InitialValue"].min()
k3_lin_range = samp_lin_max_k3 - samp_lin_min_k3

samp_lin_df = samp_lin_df.sort_values(by=['[l2]'])



#create df for normalsied rates
column_names = ["normalised lk1"]
samp_lin_rates_df = pd.DataFrame(columns = column_names)

#subtract the lower bound of the range from eahc rate in the df (normalisaiton)
samp_lin_rates_df["normalised lk1"] = samp_lin_df["Values[lk1].InitialValue"]-samp_lin_min_lk1
samp_lin_rates_df["normalised kd1"] = samp_lin_df["Values[kd1].InitialValue"]-samp_lin_min_kd1
samp_lin_rates_df["normalised k3"] = samp_lin_df["Values[k3].InitialValue"]-samp_lin_min_k3

#divide each rate cosntant column by range of rate constant (normalisation)
samp_lin_rates_df["normalised lk1"] = samp_lin_rates_df["normalised lk1"].div(lk1_lin_range)
samp_lin_rates_df["normalised kd1"] = samp_lin_rates_df["normalised kd1"].div(kd1_lin_range)
samp_lin_rates_df["normalised k3"] = samp_lin_rates_df["normalised k3"].div(k3_lin_range)

#creeate df of just species steady states
samp_lin_species_df = samp_lin_df[['[c]','[cs]','[l2]', 'Total RNA']]
#rename columns
samp_lin_species_df.columns = ['cmRNA', 'cmRNA:siRNA', 'l2RNA', 'Total RNA']

#create heatmap of normalised rates and species steady states
import seaborn as sns
fig, (ax,ax2) = plt.subplots(ncols=2, figsize=(5,10))
fig.subplots_adjust(wspace=0.02)
sns.heatmap(samp_lin_rates_df, cmap="rocket", ax=ax, cbar=False, yticklabels=False)
fig.colorbar(ax.collections[0], ax=ax,location="left", use_gridspec=False, pad=0.2)
sns.heatmap(samp_lin_species_df, cmap="Blues", ax=ax2, cbar=False, yticklabels=False)
fig.colorbar(ax2.collections[0], ax=ax2,location="right", use_gridspec=False, pad=0.2)
ax2.yaxis.tick_right()
ax.tick_params(left=False, bottom=False)
ax2.tick_params(left=False, bottom=False)
plt.show()




