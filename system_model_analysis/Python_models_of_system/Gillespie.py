#this code runs a gillespie simualtion of the reaction framework

import numpy as np
import matplotlib.pyplot as plt 


########################
####rates and species###
########################
# k1 is the rate of transcription
# kt3 is the rate of 3' processing by the 3' twister ribozyme
# kt5 is the rate of 5' processing by the 5' twister ribozyme
# k2 is the rate of circularisation
# kc is the rate of cRNA:siRNA complex formation
# kd is the rate of cRNA:siRNA complex dissocation
# k3 is the rate of dsRNAse mediated re-linearization 
# sk1 is the rate siRNA transcription 
# kd1 is the rate of linear RNA degradation
# kd2 is the rate of circular RNA degradation
# kd4 is the rate of protein degradation
# skd1 is the rate of siRNA degradation
# t1 is the rate of translation from pre-circular RNA
# t2 is the rate of translation from circular RNA
# t2 is the rate of translation from re-linearized RNA
# l1 is the number of pre-circularised linear RNAs: color = 'green', linestyle='-'
# l13 is the number of 3' processed linear RNAs: color = 'lightgreen', linestyle=':'
# l15 is the number of 3' processed linear RNAs: color = 'darkgreen', linestyle=':'
# l135 is the number of 3' & 5' processed linear RNAs: color = 'dodgerblue', linestyle=':'
# c is the number of circularised RNAs: color = 'blue', linestyle='-'
# s is the number of siRNA: color = 'red', linestyle='--'
# cs is the number of circular RNA:siRNA complexes: color = 'purple', linestyle='--'
# l2 is the number of linearised circular RNAs: color = 'cyan', linestyle='-'
# sa is the number of one half of the cleaved siRNA products: color = 'orangered', linestyle='--'
# sb is the number of other half of the cleaved siRNA products: color = 'orange', linestyle='--'
# p is the number of proteins: color = 'black', linestyle='-'
########################

#how long to run simulation for, stays outside function as refer later to calulate mean
timespan = 15000

#gillespie function
def gillespie(k1, kd1, kt3, kt5, k2, kd2, t1, t2, kd4, sk1, skd1, kc, kd, k3, t3, l1, l13, l15, l135, c, p, s, cs, l2, sa, sb, t):
    #output lists, each lsit contains species counts
    l1list = [0]
    l13list = [0]
    l15list = [0]
    l135list = [0]
    clist = [0]
    plist = [0]
    slist = [0]
    cslist =[0]
    l2list = [0]
    salist = [0]
    sblist = [0]
    tlist = [0] 
    #output is a list of each species list
    output = [tlist,l1list,l13list,l15list,l135list,clist,plist,slist,cslist,l2list,salist,sblist]
    #get species count until reach end of timespan
    while t < timespan:
           #make empty list of reaction propensities
           rlist = []
           #define reactions and what they do
           reactiona = k1 # l1+1
           reactionb = t1 * l1 # p+1
           reactionc = kd1 * l1 # l1-1
           reactiond = kt3 * l1 # l13+1, l1-1
           reactione = t1 * l13 # p+1        
           reactionf = kd1 * l13 # l13-1
           reactiong = kt5 * l1 # l15+1, l1-1
           reactionh = t1 * l15 # p+1
           reactioni = kd1 * l15 # l15-1      
           reactionj = kt5 * l13 # l135+1, l13-1
           reactionk = kt3 * l15 # l135+1, l15-1
           reactionl = t1 * l135 # p+1
           reactionm = kd1 * l135 # l135-1
           reactionn = k2 * l135 # c+1, l135-1
           reactiono = t2 * c # p+1
           reactionp = kd2 * c # c-1
           reactionq = sk1 # s+1
           reactionr = skd1 * s # s-1
           reactions = kc * c * s # cs+1, c-1, s-1
           reactiont = kd * cs # cs-1, c+1, s+1
           reactionu = k3 * cs # cs-1, l2+1, sa+1, sb+1
           reactionv = t3 * l2 # p+1
           reactionw = kd1 * l2 # l2-1
           reactionx = skd1 * sa # sa-1
           reactiony = skd1 * sb # sb-1
           reactionz = kd4 * p # p-1
           #extend empty rlist with reaction propensities
           rlist.extend(value for name, value in sorted(locals().items(), key=lambda item: item[0]) if name.startswith('reaction'))
           #calculate sum of propensities
           R = sum(rlist)
           #propensity divided by sum of propensities, so get probabiltity of reaction occuring
           xlist = [x / R for x in rlist]
           #get ranges of probabilties, so if random no. falls in a range that reaction happens
           problist = [sum(xlist[0:i+1]) for i,v in enumerate(xlist)]
           #what letter appends to what list
           l1list.append(l1)
           l13list.append(l13)
           l15list.append(l15)
           l135list.append(l135)
           clist.append(c)
           plist.append(p)
           slist.append(s)
           cslist.append(cs)
           l2list.append(l2)
           salist.append(sa)
           sblist.append(sb)
           #pick random no. between 0,1
           rand = np.random.rand()
           #if falls between x1 and x2, add/minus 1 from list/lists
           if t  > timespan*1/3: #if want to simulate introduction of siRNA
               sk1 = 0.2
           if t  > timespan*2/3: #if want to simulate introduction of siRNA
               sk1 = 0.0
           if rand >=0 and rand <= problist[0]: #reaction 1
               l1 = l1list[-1] + 1      
       	   elif rand >problist[0] and rand <= problist[1]: #reaction 2
               p = plist[-1] + 1
       	   elif rand >problist[1] and rand <= problist[2]: #reaction 3
               l1 = l1list[-1] - 1
           elif rand >problist[2] and rand <= problist[3]: #reaction 4
               l13 = l13list[-1] + 1
               l1 = l1list[-1] - 1 
       	   elif rand >problist[3] and rand <= problist[4]: #reaction 5
               p = plist[-1] + 1
           elif rand >problist[4] and rand <= problist[5]: #reaction 6
               l13 = l13list[-1] - 1 
           elif rand >problist[5] and rand <= problist[6]: #reaction 7
               l15 = l15list[-1] + 1
               l1 = l1list[-1] - 1 
           elif rand >problist[6] and rand <= problist[7]: #reaction 8
               p = plist[-1] + 1
           elif rand >problist[7] and rand <= problist[8]: #reaction 9
               l15 = l15list[-1] - 1 
           elif rand >problist[8] and rand <= problist[9]: #reaction 10
               l135 = l135list[-1] + 1
               l13 = l13list[-1] - 1 
           elif rand >problist[9] and rand <= problist[10]: #reaction 11
               l135 = l135list[-1] + 1
               l15 = l15list[-1] - 1 
           elif rand >problist[10] and rand <= problist[11]: #reaction 12
               p = plist[-1] + 1 
           elif rand >problist[11] and rand <= problist[12]: #reaction 13
               l135 = l135list[-1] - 1
           elif rand >problist[12] and rand <= problist[13]: #reaction 14
               c = clist[-1] + 1
               l135 = l135list[-1] - 1 
           elif rand >problist[13] and rand <= problist[14]: #reaction 15
               p = plist[-1] + 1 
           elif rand >problist[14] and rand <= problist[15]: #reaction 16
               c = clist[-1] - 1
           elif rand >problist[15] and rand <= problist[16]: #reaction 17
               s = slist[-1] + 1 
           elif rand >problist[16] and rand <= problist[17]: #reaction 18
               s = slist[-1] - 1 
           elif rand >problist[17] and rand <= problist[18]: #reaction 19
               cs = cslist[-1] + 1 
               c = clist[-1] - 1
               s = slist[-1] - 1
           elif rand >problist[18] and rand <= problist[19]: #reaction 20
               cs = cslist[-1] - 1 
               c = clist[-1] + 1
               s = slist[-1] + 1
           elif rand >problist[19] and rand <= problist[20]: #reaction 21
               cs = cslist[-1] - 1 
               l2 = l2list[-1] + 1
               sa = salist[-1] + 1
               sb = sblist[-1] + 1
           elif rand >problist[20] and rand <= problist[21]: #reaction 22
               p = plist[-1] + 1
           elif rand >problist[21] and rand <= problist[22]: #reaction 23
               l2 = l2list[-1] - 1
           elif rand >problist[22] and rand <= problist[23]: #reaction 24
               sa = sblist[-1] - 1
           elif rand >problist[23] and rand <= problist[24]: #reaction 25
               sb = sblist[-1] - 1
           elif rand >problist[24] and rand <= problist[25]: #reaction 26
               p = plist[-1] - 1                 
           #time increments
           tau = np.random.exponential(1/R)
           t = tlist[-1] + tau 
           tlist.append(t)
    return output


#rate constants      
k1 = 0.055
kd1 = 0.01
kt3 = 0.1
kt5 = 0.1
k2 = 0.01
kd2 = 0.001
t1 = 0.001
t2 = 0.02
kd4 = 0.01
sk1 = 0
skd1 = 0.01
kc = 0.05
kd = 0.01
k3 = 0.025
t3 = 0.001

#initial species concentrations
l1= 0
l13 = 0
l15 = 0
l135 = 0
c = 0
p = 0
t = 0
s = 0
cs = 0
l2 = 0
sa = 0
sb = 0

#no. of simulations
no_of_sim = 5
simulations = []

#run simulations
for x in range(no_of_sim):
    simulations = simulations + gillespie(k1, kd1, kt3, kt5, k2, kd2, t1, t2, kd4, sk1, skd1, kc, kd, k3, t3, l1, l13, l15, l135, c, p, s, cs, l2, sa, sb, t)
#f = number of rows generated for all simulations
f = int(len(simulations)/no_of_sim)
###########################################
#######GETTING MEAN########################
###########################################

t_range = np.arange(0, timespan)
#create an array for each list of species and time counts for each simualtion
t_array = [np.array(x) for x in simulations[::f]]
l1_array = [np.array(x) for x in simulations[1::f]]
l13_array = [np.array(x) for x in simulations[2::f]]
l15_array = [np.array(x) for x in simulations[3::f]]
l135_array = [np.array(x) for x in simulations[4::f]]
c_array = [np.array(x) for x in simulations[5::f]]
p_array = [np.array(x) for x in simulations[6::f]]
s_array = [np.array(x) for x in simulations[7::f]]
cs_array = [np.array(x) for x in simulations[8::f]]
l2_array = [np.array(x) for x in simulations[9::f]]
sa_array = [np.array(x) for x in simulations[10::f]]
sb_array = [np.array(x) for x in simulations[11::f]]

#average array of lists of species coutns fpr each simualtion
l1_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(l1_array,t_array)]) for t in t_range ]
l13_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(l13_array,t_array)]) for t in t_range ]
l15_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(l15_array,t_array)]) for t in t_range ]
l135_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(l135_array,t_array)]) for t in t_range ]
c_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(c_array,t_array)]) for t in t_range ]
p_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(p_array,t_array)]) for t in t_range ]
s_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(s_array,t_array)]) for t in t_range ]
cs_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(cs_array,t_array)]) for t in t_range ]
l2_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(l2_array,t_array)]) for t in t_range ]
sa_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(sa_array,t_array)]) for t in t_range ]
sb_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(sb_array,t_array)]) for t in t_range ]

###########################################
###################plot#####################
###########################################


counter = 0
plt.figure(figsize=(24, 8))
while counter <= ((no_of_sim*f)-1):
    linwid = 0.5
    #plot each species from each simualtion
    plt.plot(simulations[counter],simulations[counter+1], label = "Pre-circularised RNA" if counter == 0 else "", lw = linwid, color = 'green', linestyle='-') #l1
    plt.plot(simulations[counter],simulations[counter+2], label = "3'processed pre-circularised RNA" if counter == 0 else "", lw = linwid, color = 'lime', linestyle=':') #l13
    plt.plot(simulations[counter],simulations[counter+3], label = "5'processed pre-circularised RNA" if counter == 0 else "", lw = linwid, color = 'darkgreen', linestyle=':') #l15
    plt.plot(simulations[counter],simulations[counter+4], label = "3' & 5' processed pre-circularised RNA" if counter == 0 else "", lw = linwid, color = 'dodgerblue', linestyle=':') #l135
    plt.plot(simulations[counter],simulations[counter+5], label = "Circularised RNA" if counter == 0 else "", lw = linwid, color = "blue", linestyle='-') #c
    plt.plot(simulations[counter],simulations[counter+6], label = "Protein" if counter == 0 else "", lw = linwid, color = "darkgoldenrod", linestyle='-') #p
    plt.plot(simulations[counter],simulations[counter+7], label = "siRNA" if counter == 0 else "", lw = linwid, color = 'red', linestyle='--') #s
    plt.plot(simulations[counter],simulations[counter+8], label = "cRNA:siRNA complex" if counter == 0 else "", lw = linwid, color = 'purple', linestyle='--') #cs
    plt.plot(simulations[counter],simulations[counter+9], label = "Linearised RNA" if counter == 0 else "", lw = linwid, color = 'cyan', linestyle='-') #l2
    plt.plot(simulations[counter],simulations[counter+10], label = "Cleavage product a" if counter == 0 else "", lw = linwid, color = 'orange', linestyle='--') #sa
    plt.plot(simulations[counter],simulations[counter+11], label = "Cleavage product b" if counter == 0 else "", lw = linwid, color = 'orangered', linestyle='--') #sb
    #plot mean
    plt.plot(l1_average_list, lw = 1, color='black')
    plt.plot(l13_average_list, lw = 1, color='black')
    plt.plot(l15_average_list, lw = 1, color='black')
    plt.plot(l135_average_list, lw = 1, color='black')
    plt.plot(c_average_list, lw = 1, color='black')
    plt.plot(p_average_list, lw = 1, color='black') 
    plt.plot(s_average_list, lw = 1, color='black')
    plt.plot(cs_average_list, lw = 1, color='black')
    plt.plot(l2_average_list, lw = 1, color='black')
    plt.plot(sa_average_list, lw = 1, color='black')
    plt.plot(sb_average_list, lw = 1, color='black')
    plt.xlim([0, 15000])
    ax = plt.gca()
#    ax.axes.xaxis.set_visible(False)
#    ax.axes.yaxis.set_visible(False)
    plt.tick_params(axis='y', labelsize=0, length = 0)
    #ax.text(2200,9, "at t5000 siRNA induced", fontsize=15)
    #ax.text(11100,9, "at t10000 siRNA expression halted", fontsize=15)
    plt.xlabel('time',fontsize=30)
    plt.ylabel('copy number', fontsize=30)
    plt.tick_params(axis='y', labelsize=0, length = 0)
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=30)
    counter += f
    ax.set_title('Species number with time: whole system', fontsize=30)
    plt.savefig('Gillespie.jpg', dpi=300)  
plt.show()



