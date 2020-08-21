import pandas as pd

#import df and remove uneeded columns

#circualrisation scans
k1_circ_df = pd.read_csv (r'k1_circ.txt', sep='\t')
k1_circ_df = k1_circ_df.iloc[:, :-1]
k2_circ_df = pd.read_csv (r'k2_circ.txt', sep='\t')
k2_circ_df = k2_circ_df.iloc[:, :-1]
kd1_circ_df = pd.read_csv (r'kd1_circ.txt', sep='\t')
kd1_circ_df = kd1_circ_df.iloc[:, :-1]
kd2_circ_df = pd.read_csv (r'kd2_circ.txt', sep='\t')
kd2_circ_df = kd2_circ_df.iloc[:, :-1]
kt_circ_df = pd.read_csv (r'kt_circ.txt', sep='\t')
kt_circ_df = kt_circ_df.iloc[:, :-1]
kt2_circ_df = pd.read_csv (r'kt2_circ.txt', sep='\t')
kt2_circ_df = kt2_circ_df.iloc[:, :-1]

#linearisation scan
kd2_lin_df = pd.read_csv (r'kd2_lin.txt', sep='\t')
kd2_lin_df = kd2_lin_df.iloc[:, :-1]
lk1_lin_df = pd.read_csv (r'lk1_lin.txt', sep='\t')
lk1_lin_df = lk1_lin_df.iloc[:, :-1]
kd1_lin_df = pd.read_csv (r'kd1_lin.txt', sep='\t')
kd1_lin_df = kd1_lin_df.iloc[:, :-1]
k3_lin_df = pd.read_csv (r'k3_lin.txt', sep='\t')
k3_lin_df = k3_lin_df.iloc[:, :-1]


import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

###circualrisation machinery scan###
figure, axs = plt.subplots(nrows=6, ncols=6, figsize=(20, 20))
###
#colors, lineystyles & fontsizes
l1c='black'
l13c='dodgerblue'
l15c='chocolate'
l135c='goldenrod'
cc='olivedrab'
fs = 40
ktls = "-"
ktlw=1
#k1
axs[0][0].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l1]'], color=l1c)
axs[0][1].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l13]'], color=l13c)
axs[0][2].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l15]'], color=l15c)
axs[0][3].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l135]'], color=l135c)
axs[0][4].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[c]'], color=cc)
#all species
axs[0][5].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l1]'], label='l1RNA', color=l1c)
axs[0][5].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l13]'], label='l13RNA', color=l13c)
axs[0][5].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l15]'], label='l15RNA', color=l15c)
axs[0][5].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[l135]'], label='l135RNA', color=l135c)
axs[0][5].plot(k1_circ_df['# Values[k1].InitialValue'], k1_circ_df['[c]'], label='cmRNA', color=cc)
###
#kd1
axs[1][0].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l1]'], color=l1c)
axs[1][1].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l13]'], color=l13c)
axs[1][2].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l15]'], color=l15c)
axs[1][3].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l135]'], color=l135c)
axs[1][4].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[c]'], color=cc)
#all species
axs[1][5].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l1]'], color=l1c)
axs[1][5].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l13]'], color=l13c)
axs[1][5].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l15]'], color=l15c)
axs[1][5].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[l135]'], color=l135c)
axs[1][5].plot(kd1_circ_df['# Values[kd1].InitialValue'], kd1_circ_df['[c]'], color=cc)
#kt3
axs[2][0].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l1]'], color=l1c, linestyle=ktls, lw=ktlw)
axs[2][1].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l13]'], color=l13c, linestyle=ktls, lw=ktlw)
axs[2][2].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l15]'], color=l15c, linestyle=ktls, lw=ktlw)
axs[2][3].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l135]'], color=l135c, linestyle=ktls, lw=ktlw)
axs[2][4].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[c]'], color=cc, linestyle=ktls, lw=ktlw)
#all species
axs[2][5].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l1]'], color=l1c, linestyle=ktls, lw=ktlw)
axs[2][5].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l13]'], color=l13c, linestyle=ktls, lw=ktlw)
axs[2][5].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l15]'], color=l15c, linestyle=ktls, lw=ktlw)
axs[2][5].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[l135]'], color=l135c, linestyle=ktls, lw=ktlw)
axs[2][5].plot(kt_circ_df['# Values[kt3].InitialValue'], kt_circ_df['[c]'], color=cc, linestyle=ktls, lw=ktlw)
#kt5
axs[3][0].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l1]'], color=l1c, linestyle=ktls, lw=ktlw)
axs[3][1].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l13]'], color=l13c, linestyle=ktls, lw=ktlw)
axs[3][2].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l15]'], color=l15c, linestyle=ktls, lw=ktlw)
axs[3][3].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l135]'], color=l135c, linestyle=ktls, lw=ktlw)
axs[3][4].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[c]'], color=cc, linestyle=ktls, lw=ktlw)
#all species
axs[3][5].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l1]'], color=l1c, linestyle=ktls, lw=ktlw)
axs[3][5].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l13]'], color=l13c, linestyle=ktls, lw=ktlw)
axs[3][5].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l15]'], color=l15c, linestyle=ktls, lw=ktlw)
axs[3][5].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[l135]'], color=l135c, linestyle=ktls, lw=ktlw)
axs[3][5].plot(kt2_circ_df['# Values[kt5].InitialValue'], kt2_circ_df['[c]'], color=cc, linestyle=ktls, lw=ktlw)
###
#k2
axs[4][0].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l1]'], color=l1c)
axs[4][1].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l13]'], color=l13c)
axs[4][2].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l15]'], color=l15c)
axs[4][3].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l135]'], color=l135c)
axs[4][4].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[c]'], color=cc)
#all species
axs[4][5].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l1]'], color=l1c, )
axs[4][5].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l13]'], color=l13c)
axs[4][5].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l15]'], color=l15c)
axs[4][5].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[l135]'], color=l135c)
axs[4][5].plot(k2_circ_df['# Values[k2].InitialValue'], k2_circ_df['[c]'], color=cc)
###
#k2
axs[5][0].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l1]'], color=l1c)
axs[5][1].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l13]'], color=l13c)
axs[5][2].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l15]'], color=l15c)
axs[5][3].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l135]'], color=l135c)
axs[5][4].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[c]'], color=cc)
#all species
axs[5][5].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l1]'], color=l1c)
axs[5][5].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l13]'], color=l13c)
axs[5][5].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l15]'], color=l15c)
axs[5][5].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[l135]'], color=l135c)
axs[5][5].plot(kd2_circ_df['# Values[kd2].InitialValue'], kd2_circ_df['[c]'], color=cc)
###
#set column titles
axs[0][0].set_title('l1RNA', fontsize = fs)
axs[0][1].set_title('l13RNA', fontsize = fs)
axs[0][2].set_title('l15RNA', fontsize = fs)
axs[0][3].set_title('l135RNA', fontsize = fs)
axs[0][4].set_title('cmRNA', fontsize = fs)
axs[0][5].set_title('all species', fontsize = fs)
#set row titles
figure.text(-0.025, 0.90, 'k1', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.025, 0.74, 'kd1', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.025, 0.58, 'kt3', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.025, 0.42, 'kt5', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.025, 0.26, 'k2', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.025, 0.10, 'kd2', ha='center', va='center', rotation='vertical', fontsize=fs)
#change ntoaiton of y ticks
axs[1][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[1][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[1][2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[1][3].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[2][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[2][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[2][4].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[3][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[3][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[3][2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[4][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[4][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[4][2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[5][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[5][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[5][2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

#set tightness of subplots
figure.tight_layout(pad=1.25)
#set legend
figure.legend(bbox_to_anchor=(1.025, 0.55), loc='upper left', borderaxespad=0., fontsize=fs)

plt.show()


###linearisation machinery scan###
figure, axs = plt.subplots(nrows=4, ncols=4, figsize=(16, 12))
#
axs[0][2].ticklabel_format(useOffset=False, style = 'plain')
#colors, lineystyles & fontsizes
cc='green'
sc='pink'
csc='goldenrod'
l2c='dodgerblue'
fs = 40
ls = "-"
slw=1
#lk1
axs[0][0].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[c]'], color=cc, label='cmRNA')
axs[0][1].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[cs]'], color=csc, label = 'cmRNA:siRNA')
axs[0][1].ticklabel_format(useOffset=False, style = 'plain')
axs[0][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[0][2].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[l2]'], color=l2c, label = 'l2rna')
#all species
axs[0][3].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[c]'], color=cc)
axs[0][3].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[cs]'], color=csc)
axs[0][3].plot(lk1_lin_df['# Values[lk1].InitialValue'], lk1_lin_df['[l2]'], color=l2c)
###
#kd1
axs[1][0].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[c]'], color=cc)
axs[1][1].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[cs]'], color=csc)
axs[1][1].ticklabel_format(useOffset=False, style = 'plain')
axs[1][2].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[l2]'], color=l2c)
#all species
axs[1][3].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[c]'], color=cc)
axs[1][3].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[cs]'], color=csc)
axs[1][3].plot(kd1_lin_df['# Values[lkd1].InitialValue'], kd1_lin_df['[l2]'], color=l2c)
#k3
axs[2][0].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[c]'], color=cc)
axs[2][1].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[cs]'], color=csc)
axs[2][2].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[l2]'], color=l2c)
#all species
axs[2][3].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[c]'], color=cc)
axs[2][3].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[cs]'], color=csc)
axs[2][3].plot(k3_lin_df['# Values[k3].InitialValue'], k3_lin_df['[l2]'], color=l2c)
#kd1
axs[3][0].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[c]'], color=cc)
axs[3][0].ticklabel_format(useOffset=False, style = 'plain')
axs[3][1].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[cs]'], color=csc)
axs[3][2].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[l2]'], color=l2c)
axs[3][2].ticklabel_format(useOffset=False, style = 'plain')

#all species
axs[3][3].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[c]'], color=cc)
axs[3][3].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[cs]'], color=csc)
axs[3][3].plot(kd2_lin_df['# Values[kd2].InitialValue'], kd2_lin_df['[l2]'], color=l2c)
#k3
###
#set column titles
axs[0][0].set_title('cmRNA', fontsize = fs, y=1.08)
axs[0][1].set_title('cmRNA:siRNA', fontsize = fs, y=1.08)
axs[0][2].set_title('l2RNA', fontsize = fs, y=1.08)
axs[0][3].set_title('all species', fontsize = fs, y=1.08)
# #set row titles
figure.text(-0.02, 0.80, 'sk1', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.02, 0.62, 'kd1', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.02, 0.38, 'k3', ha='center', va='center', rotation='vertical', fontsize=fs)
figure.text(-0.02, 0.15, 'kd2', ha='center', va='center', rotation='vertical', fontsize=fs)
#change ntoaiton of y ticks
axs[0][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[1][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[0][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[1][1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[2][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[3][0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
axs[3][2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#set tightness of subplots
figure.tight_layout(pad=1.1)
#set legend
figure.legend(bbox_to_anchor=(1.05, 0.6), loc='upper left', borderaxespad=0., fontsize=fs)

plt.show()


