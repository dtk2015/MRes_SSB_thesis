import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

#this script correlats ss% and mfe from the kim 2008 study agaisnt half life

five = pd.read_csv("5UTR_nupack.csv")
#drop rows where no secondary structure was predicted
five = five.dropna(axis=0, subset=['mfe_secondary_structure'])
#drop rows if length < 10
five = five[five.length > 10]
#function to get gc%
def gc_perc(seq):
    count = seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c') 
    return round((count / len(seq))*100, 2)

#add gc% as a column to each df
five['gc%'] = five['sequence'].apply(gc_perc)

def sec_str(str):
    count = str.count('(') + str.count(')')
    return round((count / len(str))*100, 2)

#add sec_str% as a column to each df
five['sec_str%'] = five['mfe_secondary_structure'].apply(sec_str)


five.to_csv("5mfe_nupack.csv", mode = 'w', index=False)

##############################################################################
#########################half lives & cds#####################################
##############################################################################

hl = pd.read_csv("half_life.csv")

hl.columns = ['gene_name',
 'gene',
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

five = five.merge(hl,on='gene_name',how='left') 

five = five.dropna(how='all', subset=['mRNA half-lifea in min at µ=0.10 h-1',
 'mRNA half-lifea in min at µ=0.20 h-1',
 'mRNA half-lifea in min at µ=0.40 h-1',
 'mRNA half-lifea in min at µ=0.63 h-1']) 


##############################################################################
########################### 5'plots & stats ##################################
##############################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.stats as stats

#mfe
y5mfe = np.array(five['mfe'])
x51 = np.array(five['mRNA half-lifea in min at µ=0.10 h-1']).reshape((-1, 1))
model5mfe1 = LinearRegression()
model5mfe1.fit(x51, y5mfe)
model5mfe1 = LinearRegression().fit(x51, y5mfe)
r_sqy5mfex51 = model5mfe1.score(x51, y5mfe)
print('coefficient of determination for 5 UTR (mfe) at u=0.1:', r_sqy5mfex51)


x52 = np.array(five['mRNA half-lifea in min at µ=0.20 h-1']).reshape((-1, 1))
model5mfe2 = LinearRegression()
model5mfe2.fit(x52, y5mfe)
model5mfe2 = LinearRegression().fit(x52, y5mfe)
r_sqy5mfex52 = model5mfe2.score(x52, y5mfe, )
print('coefficient of determination for 5 UTR (mfe) at  u=0.2:', r_sqy5mfex52)


x54 = np.array(five['mRNA half-lifea in min at µ=0.40 h-1']).reshape((-1, 1))
model5mfe4 = LinearRegression()
model5mfe4.fit(x54, y5mfe)
model5mfe4 = LinearRegression().fit(x54, y5mfe)
r_sqy5mfex54 = model5mfe4.score(x54, y5mfe)
print('coefficient of determination for 5 UTR (mfe) at  u=0.4:', r_sqy5mfex54)


x563 = np.array(five['mRNA half-lifea in min at µ=0.63 h-1']).reshape((-1, 1))
model5mfe63 = LinearRegression()
model5mfe63.fit(x563, y5mfe)
model5mfe63 = LinearRegression().fit(x563, y5mfe)
r_sqy5mfex563 = model5mfe63.score(x563, y5mfe,)
print('coefficient of determination for 5 UTR (mfe) at  u=0.63:', r_sqy5mfex563)


#ss%
y5ss = np.array(five['sec_str%'])
model5ss1 = LinearRegression()
model5ss1.fit(x51, y5ss)
model5ss1 = LinearRegression().fit(x51, y5ss)
model5ss1 = LinearRegression().fit(x563, y5ss)
r_sqy5ssx51 = model5ss1.score(x51, y5ss)
print('coefficient of determination for 5 UTR (ss%) at  u=0.10:', r_sqy5ssx51)

model5ss2 = LinearRegression()
model5ss2.fit(x52, y5ss)
model5ss2 = LinearRegression().fit(x52, y5ss)
r_sqy5ssx52 = model5ss2.score(x52, y5ss)
print('coefficient of determination for 5 UTR (ss%) at  u=0.20:', r_sqy5ssx52)

model5ss4 = LinearRegression()
model5ss4.fit(x54, y5ss)
model5ss4 = LinearRegression().fit(x54, y5ss)
r_sqy5ssx54 = model5ss4.score(x54, y5ss)
print('coefficient of determination for 5 UTR (ss%) at  u=0.40:', r_sqy5ssx54)

model5ss63 = LinearRegression()
model5ss63.fit(x563, y5ss)
model5ss63 = LinearRegression().fit(x563, y5ss)
r_sqy5ssx563 = model5ss63.score(x563, y5ss)
print('coefficient of determination for 5 UTR (ss%) at  u=0.63:', r_sqy5ssx563)

#ss%
y5ss = np.array(five['sec_str%'])
model5ss1 = LinearRegression()
model5ss1.fit(x51, y5ss)
model5ss1 = LinearRegression().fit(x51, y5ss)

model5ss2 = LinearRegression()
model5ss2.fit(x52, y5ss)
model5ss2 = LinearRegression().fit(x52, y5ss)

model5ss4 = LinearRegression()
model5ss4.fit(x54, y5ss)
model5ss4 = LinearRegression().fit(x54, y5ss)

model5ss63 = LinearRegression()
model5ss63.fit(x563, y5ss)
model5ss63 = LinearRegression().fit(x563, y5ss)

print(r_sqy5mfex51)
print(r_sqy5mfex52)
print(r_sqy5mfex54)
print(r_sqy5ssx563)
print(r_sqy5ssx51)
print(r_sqy5ssx52)
print(r_sqy5ssx54)
print(r_sqy5ssx563)
print(len(x51))
print(len(x52))
print(len(x54))
print(len(x563))
print(len(y5mfe))
print(len(y5ss))

#
ls = "-"
lw = 1
rc = "r"
spc = "dodgerblue"
sps = 1
fs = "15"
#
figure, axs = plt.subplots(nrows=2, ncols=4, figsize=(10, 5))
ylab0 = 'mfe'
ylab1 = 'secondary structure %'
xlab = 'half life (min)'

axs[0][0].set_title('µ=0.10 h-1', fontsize = fs)
axs[0][1].set_title('µ=0.20 h-1', fontsize = fs)
axs[0][2].set_title('µ=0.40 h-1', fontsize = fs)
axs[0][3].set_title('µ=0.63 h-1', fontsize = fs)

axs[0][0].set_xlabel(xlab)
axs[0][1].set_xlabel(xlab)
axs[0][2].set_xlabel(xlab)
axs[0][3].set_xlabel(xlab)
axs[1][0].set_xlabel(xlab)
axs[1][1].set_xlabel(xlab)
axs[1][2].set_xlabel(xlab)
axs[1][3].set_xlabel(xlab)



figure.text(-0.0, 0.75, 'mfe', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.0, 0.3, 'secondary structure %', ha='center', va='center', rotation='vertical', fontsize=fs)

axs[0][0].scatter(x51, y5mfe,s=sps)
axs[0][0].plot(x51, model5mfe1.predict(x51),color=rc, lw = lw, linestyle = ls)
axs[0][1].scatter(x52, y5mfe,s=sps)
axs[0][1].plot(x52, model5mfe2.predict(x52),color=rc, lw = lw, linestyle = ls)
axs[0][2].scatter(x54, y5mfe,s=sps)
axs[0][2].plot(x54, model5mfe4.predict(x54),color=rc, lw = lw, linestyle = ls)
axs[0][3].scatter(x563, y5mfe,s=sps)
axs[0][3].plot(x563, model5mfe63.predict(x563),color=rc, lw = lw, linestyle = ls)

axs[1][0].scatter(x51, y5ss,s=sps)
axs[1][0].plot(x51, model5ss1.predict(x51),color=rc, lw = lw, linestyle = ls)
axs[1][1].scatter(x52, y5ss,s=sps)
axs[1][1].plot(x52, model5ss2.predict(x52),color=rc, lw = lw, linestyle = ls)
axs[1][2].scatter(x54, y5ss,s=sps)
axs[1][2].plot(x54, model5ss4.predict(x54),color=rc, lw = lw, linestyle = ls)
axs[1][3].scatter(x563, y5ss,s=sps)
axs[1][3].plot(x563, model5ss63.predict(x563),color=rc, lw = lw, linestyle = ls)

figure.tight_layout(pad=1.25)
plt.show()





x51 = np.array(five['mRNA half-lifea in min at µ=0.10 h-1'])
x52 = np.array(five['mRNA half-lifea in min at µ=0.20 h-1'])
x54 = np.array(five['mRNA half-lifea in min at µ=0.40 h-1'])
x563 = np.array(five['mRNA half-lifea in min at µ=0.63 h-1'])
y5mfe = np.array(five['mfe'])
y5ss = np.array(five['sec_str%'])



print("5' UTR mfe vs half life at µ=0.10 h-1: cod & p-value of cod:", stats.pearsonr(x51, y5mfe))
print("5' UTR secondary structure % vs half life at µ=0.10 h-1: cod & p-value of cod:", stats.pearsonr(x51, y5ss))
print("5' UTR mfe vs half life at µ=0.20 h-1: cod & p-value of cod:", stats.pearsonr(x52, y5mfe))
print("5' UTR secondary structure % vs half life at µ=0.10 h-1: cod & p-value of cod:", stats.pearsonr(x52, y5ss))
print("5' UTR mfe vs half life at µ=0.40 h-1: cod & p-value of cod:", stats.pearsonr(x54, y5mfe))
print("5' UTR secondary structure % vs half life at µ=0.10 h-1: cod & p-value of cod:", stats.pearsonr(x54, y5ss))
print("5' UTR mfe vs half life at µ=0.63 h-1: cod & p-value of cod:", stats.pearsonr(x563, y5mfe))
print("5' UTR secondary structure % vs half life at µ=0.10 h-1: cod & p-value of cod:", stats.pearsonr(x563, y5ss))

print("")
print(stats.pearsonr(x51, y5mfe))
print(stats.pearsonr(x52, y5mfe))
print(stats.pearsonr(x54, y5mfe))
print(stats.pearsonr(x563, y5mfe))
print(stats.pearsonr(x51, y5ss))
print(stats.pearsonr(x52, y5ss))
print(stats.pearsonr(x54, y5ss))
print(stats.pearsonr(x563, y5ss))


