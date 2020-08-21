import pandas as pd
import re
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

#this code correlates TIR with half life 

######################################################################
######################################################################
######################################################################
#open df containing all traits including TIR
full_df = pd.read_csv("RBS_calc.csv")

##############################################################################
#########################half lives & cds#####################################
###############################################################################

#load half lfie data, rename columsn and drop unnecessary rows
hl = pd.read_csv("half_life_esquerre.csv")

hl.columns = ['name',
 'gene_name',
 'mRNA concentrationa in AU/µg DCW at µ=0.10 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.20 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.40 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.63 h-1',
 'mRNA half-lifea in min at µ=0.10 h-1',
 'mRNA half-lifea in min at µ=0.20 h-1',
 'mRNA half-lifea in min at µ=0.40 h-1',
 'mRNA half-lifea in min at µ=0.63 h-1',
 'Gene in Operonb?',
 'Number of genes in operonb',
 'Position in operonb',
 'Gene is essentialc? ',
 'COGd',
 'Cell location of gene productd',
 'length of ORFe in bp',
 '%GC of ORF',
 "length of 5'UTRf in bp",
 "%GC of 5'UTR",
 "length of ORF+5'UTR in bp",
 "%GC of ORF+5'UTR",
 'CAI of ORF',
 'Protein with signal peptided?']


hl = hl.drop([
 'mRNA concentrationa in AU/µg DCW at µ=0.10 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.20 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.40 h-1',
 'mRNA concentrationa in AU/µg DCW at µ=0.63 h-1',
 'Gene in Operonb?',
 'Number of genes in operonb',
 'Position in operonb',
 'Gene is essentialc? ',
 'COGd',
 'Cell location of gene productd',
 'length of ORFe in bp',
 '%GC of ORF',
 "length of 5'UTRf in bp",
 "%GC of 5'UTR",
 "length of ORF+5'UTR in bp",
 "%GC of ORF+5'UTR",
 'CAI of ORF',
 'Protein with signal peptided?']
    , axis=1)

#merge HL data with RBS calc data
df_full = full_df.merge(hl,on='name',how='left') 

#drop nan values
df_full = df_full.dropna(how='all', subset=['mRNA half-lifea in min at µ=0.10 h-1',
 'mRNA half-lifea in min at µ=0.20 h-1',
 'mRNA half-lifea in min at µ=0.40 h-1',
 'mRNA half-lifea in min at µ=0.63 h-1']) 

df_full = df_full[df_full['TIR'].notna()]

##############################################################################
######################### plots & stats ######################################
##############################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.stats as stats

#calauclte lienar regressison and coefficient of determiantion
y = np.array(df_full['TIR'])
#calcualting lr & cod for eahc condition
x1 = np.array(df_full['mRNA half-lifea in min at µ=0.10 h-1']).reshape((-1, 1))
model1 = LinearRegression()
model1.fit(x1, y)
model1 = LinearRegression().fit(x1, y)
r_sqyx1 = model1.score(x1, y)
print('coefficient of determination of u=0.1:', r_sqyx1)

#calcualting lr & cod for eahc condition
x2 = np.array(df_full['mRNA half-lifea in min at µ=0.20 h-1']).reshape((-1, 1))
model2 = LinearRegression()
model2.fit(x2, y)
model2 = LinearRegression().fit(x2, y)
r_sqyx2 = model2.score(x2, y)
print('coefficient of determination of u=0.2:', r_sqyx2)

#calcualting lr & cod for eahc condition
x4 = np.array(df_full['mRNA half-lifea in min at µ=0.40 h-1']).reshape((-1, 1))
model4 = LinearRegression()
model4.fit(x4, y)
model4 = LinearRegression().fit(x4, y)
r_sqyx4 = model4.score(x4, y)
print('coefficient of determination of u=0.4:', r_sqyx4)

#calcualting lr & cod for eahc condition
x63 = np.array(df_full['mRNA half-lifea in min at µ=0.63 h-1']).reshape((-1, 1))
model63 = LinearRegression()
model63.fit(x63, y)
model63 = LinearRegression().fit(x63, y)
r_sqyx63 = model63.score(x63, y)
print('coefficient of determination of u=0.63:', r_sqyx63)
#
#customisign aspects of plot
ls = ":"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 0.1
fs = "20"
###########
###plots###
###########
figure, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 5))
ylab = 'TIR'
xlab = 'half life (min)'
#set title
axs[0][0].set_title('TIR vs half life at µ=0.10 h-1')
axs[0][1].set_title('TIR vs half life at µ=0.20 h-1')
axs[1][0].set_title('TIR vs half life at µ=0.40 h-1')
axs[1][1].set_title('TIR vs half life at µ=0.63 h-1')
#set label axis'
axs[0][0].set_xlabel(xlab)
axs[1][0].set_xlabel(xlab)
axs[0][1].set_xlabel(xlab)
axs[1][1].set_xlabel(xlab)
axs[0][0].set_ylabel(ylab)
axs[0][1].set_ylabel(ylab)
axs[1][0].set_ylabel(ylab)
axs[1][1].set_ylabel(ylab)


#plot data poitns and lr
axs[0][0].scatter(x1, y,s=sps)
axs[0][0].plot(x1, model1.predict(x1),color=rc, lw = lw, linestyle = ls)
axs[0][1].scatter(x2, y,s=sps)
axs[0][1].plot(x2, model2.predict(x2),color=rc, lw = lw, linestyle = ls)
axs[1][0].scatter(x4, y,s=sps)
axs[1][0].plot(x4, model4.predict(x4),color=rc, lw = lw, linestyle = ls)
axs[1][1].scatter(x63, y,s=sps)
axs[1][1].plot(x63, model63.predict(x63),color=rc, lw = lw, linestyle = ls)

figure.tight_layout(pad=1.25)
plt.show()


#calcualte coeffcieint of detemriantion & p-value
x1 = np.array(df_full['mRNA half-lifea in min at µ=0.10 h-1'])
x2 = np.array(df_full['mRNA half-lifea in min at µ=0.20 h-1'])
x4 = np.array(df_full['mRNA half-lifea in min at µ=0.40 h-1'])
x63 = np.array(df_full['mRNA half-lifea in min at µ=0.63 h-1'])
print('TIR vs half life at µ=0.10 h-1: cod & p-value of cod:', stats.pearsonr(x1, y))
print('TIR vs half life at µ=0.20 h-1: cod & p-value of cod:', stats.pearsonr(x2, y))
print('TIR vs half life at µ=0.40 h-1: cod & p-value of cod:', stats.pearsonr(x4, y))
print('TIR vs half life at µ=0.63 h-1: cod & p-value of cod:', stats.pearsonr(x63, y))

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

#make list of TIR outliers
outliers_TIR = detect_outlier(df_full['TIR'])
#create a df of these outliers and mark eac as treue
outliers_TIR_df =pd.DataFrame(list(zip(outliers_TIR, [True]*len(outliers_TIR))), columns=["TIR", "outlier_TIR"])
#malke list of 
outliers_hl_10 = detect_outlier(df_full['mRNA half-lifea in min at µ=0.10 h-1'])
outliers_hl_10_df =pd.DataFrame(list(zip(outliers_hl_10, [True]*len(outliers_hl_10))), columns=["mRNA half-lifea in min at µ=0.10 h-1", "outlier_mRNA half-lifea in min at µ=0.10 h-1"])
outliers_hl_20 = detect_outlier(df_full['mRNA half-lifea in min at µ=0.20 h-1'])
outliers_hl_20_df =pd.DataFrame(list(zip(outliers_hl_20, [True]*len(outliers_hl_20))), columns=["mRNA half-lifea in min at µ=0.20 h-1", "outlier_mRNA half-lifea in min at µ=0.20 h-1"])
outliers_hl_40 = detect_outlier(df_full['mRNA half-lifea in min at µ=0.40 h-1'])
outliers_hl_40_df =pd.DataFrame(list(zip(outliers_hl_40, [True]*len(outliers_hl_40))), columns=["mRNA half-lifea in min at µ=0.40 h-1", "outlier_mRNA half-lifea in min at µ=0.40 h-1"])
outliers_hl_63 = detect_outlier(df_full['mRNA half-lifea in min at µ=0.63 h-1'])
outliers_hl_63_df =pd.DataFrame(list(zip(outliers_hl_63, [True]*len(outliers_hl_63))), columns=["mRNA half-lifea in min at µ=0.63 h-1", "outlier_mRNA half-lifea in min at µ=0.63 h-1"])



#create individual df of just half lives and TIR
hl_10 = df_full[['mRNA half-lifea in min at µ=0.10 h-1']]
hl_20 = df_full[['mRNA half-lifea in min at µ=0.20 h-1']]
hl_40 = df_full[['mRNA half-lifea in min at µ=0.40 h-1']]
hl_63 = df_full[['mRNA half-lifea in min at µ=0.63 h-1']]
tir = df_full[['TIR']]

#merge with outlier data
hl_10 = hl_10.merge(outliers_hl_10_df,on='mRNA half-lifea in min at µ=0.10 h-1',how='left') 
hl_20 = hl_20.merge(outliers_hl_20_df,on='mRNA half-lifea in min at µ=0.20 h-1',how='left') 
hl_40 = hl_40.merge(outliers_hl_40_df,on='mRNA half-lifea in min at µ=0.40 h-1',how='left') 
hl_63 = hl_63.merge(outliers_hl_63_df,on='mRNA half-lifea in min at µ=0.63 h-1',how='left') 
tir = tir.merge(outliers_TIR_df,on='TIR',how='left') 

#drop rows that are outliers
hl_10.drop(hl_10[hl_10['outlier_mRNA half-lifea in min at µ=0.10 h-1'] == True].index, inplace = True) 
hl_20.drop(hl_20[hl_20['outlier_mRNA half-lifea in min at µ=0.20 h-1'] == True].index, inplace = True) 
hl_40.drop(hl_40[hl_40['outlier_mRNA half-lifea in min at µ=0.40 h-1'] == True].index, inplace = True) 
hl_63.drop(hl_63[hl_63['outlier_mRNA half-lifea in min at µ=0.63 h-1'] == True].index, inplace = True) 
tir.drop(tir[tir['outlier_TIR'] == True].index, inplace = True) 

#merge with full dataframe
hl_10 = df_full[df_full['mRNA half-lifea in min at µ=0.10 h-1'].isin(hl_10['mRNA half-lifea in min at µ=0.10 h-1'])]
hl_20 = df_full[df_full['mRNA half-lifea in min at µ=0.20 h-1'].isin(hl_20['mRNA half-lifea in min at µ=0.20 h-1'])]
hl_40 = df_full[df_full['mRNA half-lifea in min at µ=0.40 h-1'].isin(hl_40['mRNA half-lifea in min at µ=0.40 h-1'])]
hl_63 = df_full[df_full['mRNA half-lifea in min at µ=0.63 h-1'].isin(hl_63['mRNA half-lifea in min at µ=0.63 h-1'])]
tir = df_full[df_full['TIR'].isin(tir['TIR'])]

#subset the dataframe into just the no outlier value and name 
hl_10_no = hl_10[["name", "TIR", "mRNA half-lifea in min at µ=0.10 h-1"]]
hl_20_no = hl_20[["name", "TIR", "mRNA half-lifea in min at µ=0.20 h-1"]]
hl_40_no = hl_40[["name", "TIR", "mRNA half-lifea in min at µ=0.40 h-1"]]
hl_63_no = hl_63[["name", "TIR", "mRNA half-lifea in min at µ=0.63 h-1"]]


##############################################################################
######################### plots & stats ######################################
##############################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.stats as stats

#calcualtign lr and coefficent of detemriantion
y1_no = np.array(hl_10_no['TIR'])
x1_no = np.array(hl_10_no['mRNA half-lifea in min at µ=0.10 h-1']).reshape((-1, 1))
model1no = LinearRegression()
model1no.fit(x1_no, y1_no)
model1no = LinearRegression().fit(x1_no, y1_no)
r_sqy1x1_no = model1no.score(x1_no, y1_no)
print('coefficient of determination of u=0.1:', r_sqy1x1_no)

y2_no = np.array(hl_20_no['TIR'])
x2_no = np.array(hl_20_no['mRNA half-lifea in min at µ=0.20 h-1']).reshape((-1, 1))
model2no = LinearRegression()
model2no.fit(x2_no, y2_no)
model2no = LinearRegression().fit(x2_no, y2_no)
r_sqy2x2_no = model2no.score(x2_no, y2_no)
print('coefficient of determination of u=0.2:', r_sqy2x2_no)

y4_no = np.array(hl_40_no['TIR'])
x4_no = np.array(hl_40_no['mRNA half-lifea in min at µ=0.40 h-1']).reshape((-1, 1))
model4no = LinearRegression()
model4no.fit(x4_no, y4_no)
model4no = LinearRegression().fit(x4_no, y4_no)
r_sqy4x4_no = model4no.score(x4_no, y4_no)
print('coefficient of determination of u=0.4:', r_sqy4x4_no)


y63_no = np.array(hl_63_no['TIR'])
x63_no = np.array(hl_63_no['mRNA half-lifea in min at µ=0.63 h-1']).reshape((-1, 1))
model63no = LinearRegression()
model63no.fit(x63_no, y63_no)
model63no = LinearRegression().fit(x63_no, y63_no)
r_sqy63x63_no = model63no.score(x63_no, y63_no)
print('coefficient of determination of u=0.63:', r_sqy63x63_no)

#customisign aspects of plot
ls = ":"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 0.1
fs = "20"
###########
###plots###
###########
ls = ":"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 0.1
fs = "20"
#
figure, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 5))
ylab = 'TIR'
xlab = 'half life (min)'
#set title
axs[0][0].set_title('TIR vs half life at µ=0.10 h-1')
axs[0][1].set_title('TIR vs half life at µ=0.20 h-1')
axs[1][0].set_title('TIR vs half life at µ=0.40 h-1')
axs[1][1].set_title('TIR vs half life at µ=0.63 h-1')
#set axis' titles
axs[0][0].set_xlabel(xlab)
axs[1][0].set_xlabel(xlab)
axs[0][1].set_xlabel(xlab)
axs[1][1].set_xlabel(xlab)
axs[0][0].set_ylabel(ylab)
axs[0][1].set_ylabel(ylab)
axs[1][0].set_ylabel(ylab)
axs[1][1].set_ylabel(ylab)

#plot datapoitns and lr
axs[0][0].scatter(x1_no, y1_no,s=sps)
axs[0][0].plot(x1_no, model1no.predict(x1_no),color=rc, lw = lw, linestyle = ls)
axs[0][1].scatter(x2_no, y2_no,s=sps)
axs[0][1].plot(x2_no, model2no.predict(x2_no),color=rc, lw = lw, linestyle = ls)
axs[1][0].scatter(x4_no, y4_no,s=sps)
axs[1][0].plot(x4_no, model4no.predict(x4_no),color=rc, lw = lw, linestyle = ls)
axs[1][1].scatter(x63_no, y63_no,s=sps)
axs[1][1].plot(x63_no, model63no.predict(x63_no),color=rc, lw = lw, linestyle = ls)

figure.tight_layout(pad=1.25)
plt.show()

#calcualte coeffcieint of detemriantion & p-value
x1_no = np.array(hl_10_no['mRNA half-lifea in min at µ=0.10 h-1'])
x2_no = np.array(hl_20_no['mRNA half-lifea in min at µ=0.20 h-1'])
x4_no = np.array(hl_40_no['mRNA half-lifea in min at µ=0.40 h-1'])
x63_no = np.array(hl_63_no['mRNA half-lifea in min at µ=0.63 h-1'])
y1_no = np.array(hl_10_no['TIR'])
y2_no = np.array(hl_20_no['TIR'])
y4_no = np.array(hl_40_no['TIR'])
y63_no = np.array(hl_63_no['TIR'])

print('TIR vs half life at µ=0.10 h-1: cod & p-value of cod:', stats.pearsonr(x1_no, y1_no))
print('TIR vs half life at µ=0.20 h-1: cod & p-value of cod:', stats.pearsonr(x2_no, y2_no))
print('TIR vs half life at µ=0.40 h-1: cod & p-value of cod:', stats.pearsonr(x4_no, y4_no))
print('TIR vs half life at µ=0.63 h-1: cod & p-value of cod:', stats.pearsonr(x63_no, y63_no))
