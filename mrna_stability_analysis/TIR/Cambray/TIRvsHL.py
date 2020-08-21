import pandas as pd
import re
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

#this code plots half life agaisnt TIR from Cambrays 244,000 sytnehtic trnascipt library

######################################################################
######################################################################
######################################################################
#open df containing all traits including TIR
rbs_calc = pd.read_csv("RBS_calc.csv")

##############################################################################
#########################half lives & cds#####################################
##############################################################################

hl = pd.read_csv("pheno_meas.csv")

#merge
df_full = rbs_calc.merge(hl,on='name',how='left') 

filtered_df = df_full[['name', 'tir', 'halflife.rna.dna.mean']]
filtered_df = filtered_df.dropna(axis='rows')

##############################################################################
######################### plots & stats ######################################
##############################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.stats as stats

y = np.array(filtered_df['tir'])
#calcualting lr & cod for eahc condition
x = np.array(filtered_df['halflife.rna.dna.mean']).reshape((-1, 1))
model = LinearRegression()
model.fit(x, y)
model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print('coefficient of determination:', r_sq)
#
ls = ":"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 0.1
fs = "20"
#
a = plt.title('TIR vs half life')
a = plt.ylabel("TIR")
a = plt.xlabel("half life (min)")
a = plt.scatter(x, y,s=sps)
a = plt.plot(x, model.predict(x),color=rc, lw = lw, linestyle = ls)
a = plt.ticklabel_format(useOffset=False, style = 'plain')
plt.show()

#calcualte coeffcieint of detemriantion & p-value
x1 = np.array(filtered_df['halflife.rna.dna.mean'])
y1 = np.array(filtered_df['tir'])

print('TIR vs half life: cod & p-value of cod:', stats.pearsonr(x1, y1))


##############################################################################
####################### removed outliers #####################################
##############################################################################
print("")
print("#################################################")
print("############### removing outliers ###############")
print("#################################################")
print("")
#calcualte outliers based on z score
def detect_outlier(data):    
    outliers= []
    threshold=3
    mean = np.mean(data)
    std =np.std(data)
        
    for y in data:
        z_score= (y - mean)/std 
        if np.abs(z_score) > threshold:
            outliers.append(y)
    return outliers

#not removing outliers for TIR
tir = filtered_df[['name', 'tir']]




outliers_hl = detect_outlier(filtered_df['halflife.rna.dna.mean'])
#create a df of these outliers and mark eac as treue
outliers_hl_df =pd.DataFrame(list(zip(outliers_hl, [True]*len(outliers_hl))), columns=["halflife.rna.dna.mean", "outlier_halflife.rna.dna.mean"])
#malke list of 


#create individual df of just half lives and TIR
hl = filtered_df[['name', 'halflife.rna.dna.mean']]

#merge with outlier data
hl = hl.merge(outliers_hl_df,on='halflife.rna.dna.mean',how='left') 

#drop rows that are outliers
hl.drop(hl[hl['outlier_halflife.rna.dna.mean'] == True].index, inplace = True) 

#
#merge with outlier data
hl_no = hl.merge(tir,on='name',how='left') 

#drop rows with no tir
hl_no = hl_no.dropna(subset = ['tir'])

##############################################################################
######################### plots & stats ######################################
##############################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.stats as stats

y1 = np.array(hl_no['tir'])
#calcualting lr & cod for eahc condition
x1 = np.array(hl_no['halflife.rna.dna.mean']).reshape((-1, 1))
modelno = LinearRegression()
modelno.fit(x1, y1)
modelno = LinearRegression().fit(x1, y1)
r_sq1 = model.score(x1, y1)
print('coefficient of determination:', r_sq1)
#
ls = ":"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 0.1
fs = "20"
#
b = plt.title('TIR vs half life')
b = plt.ylabel("TIR")
b = plt.xlabel("half life (min)")
b = plt.scatter(x1,y1,s=sps)
b = plt.plot(x1, modelno.predict(x1),color=rc, lw = lw, linestyle = ls)
b = plt.ticklabel_format(useOffset=False, style = 'plain')
plt.show()

#calcualte coeffcieint of detemriantion & p-value
x1_no = np.array(hl_no['halflife.rna.dna.mean'])
y1_no = np.array(hl_no['tir'])

print('TIR vs half life: cod & p-value of cod:', stats.pearsonr(x1_no, y1_no))

#

figure, axs = plt.subplots(2, figsize=(10, 5))
axs[0].set_title('TIR vs half life')
axs[1].set_title('TIR vs half life (outliers removed)')
axs[0].set_xlabel("half life (mins)")
axs[1].set_xlabel("half life (mins")
axs[0].set_ylabel("TIR")
axs[1].set_ylabel("TIR")
axs[0].scatter(x, y,s=sps)
axs[0].plot(x, model.predict(x),color=rc, lw = lw, linestyle = ls)
axs[1].scatter(x1, y1,s=sps)
axs[1].plot(x1, modelno.predict(x1),color=rc, lw = lw, linestyle = ls)
axs[0].ticklabel_format(useOffset=False, style = 'plain')
axs[1].ticklabel_format(useOffset=False, style = 'plain')
figure.tight_layout(pad=1.25)
plt.show()


