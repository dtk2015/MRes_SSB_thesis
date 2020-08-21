#this code runs a gillespie simualtion of the reaction framework

import numpy as np
import scipy.integrate
from timeit import default_timer as timer
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

# set rate constants
k1 = 0.055
kd1 = 0.01
kt3 = 0.1
kt5 = 0.1
k2 = 0.01
kd2 = 0.001
t1 = 0.001
t2 = 0.02
kd4 = 0.01
sk = 0
skd = 0.01
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
###########################################
###########Whole model#####################
###########################################

#solving ODEs, and obtainign determisntic soltuion for eahc species with time
def drdt(t,r,k1,kt3,kt5,k2,kc,kd,k3,sk,skd,kd1,kd2,kd4,t1,t2,t3):
    l1, l13, l15, l135, c, s, cs, l2, sa, sb, p = r
    #specify at what time siRNA expression is induced
    if t > 5000:
        sk = 0.2
    #specify at what time siRNA expression is halted
    if t > 10000:
        sk = 0.0
    #ODEs
    dl1dt = k1-((kt5+kt3+kd1)*l1)
    dl13dt = kt3*l1 - kt5*l13
    dl15dt = kt5*l1 - kt3*l15
    dl135dt = ((kt5*l13)+(kt3*l15))-((k2+kd1)*l135)
    dcdt = ((k2*l135)+(kd*cs)) - ((kd2*c)+(kc*(c*s)))
    dsdt = (sk+(kd*cs))-((kd1*s) + (kc*(c*s)))
    dcsdt = kc*(c*s) - ((cs*k3)+(cs*kd))
    dl2dt = (cs*kc) - (l2*kd1)
    dsadt = (cs*kc) - (sa*skd)
    dsbdt = (cs*kc) - (sb*skd)
    dpdt = ((t1*(l1+l13+l15+l135)) + (t2*c) + (t3*l2))-(kd4*p)
    return (dl1dt, dl13dt, dl15dt, dl135dt, dcdt, dsdt, dcsdt, dl2dt, dsadt, dsbdt, dpdt)


#initial species count
l10, l130, l150, l1350, c0, s0, cs0, l20, sa0, sb0, p0 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
drdt_withks = lambda t,r: drdt(t,r,k1,kt3,kt5,k2,kc,kd,k3,sk,skd,kd1,kd2,kd4,t1,t2,t3)

#solve odes
start = timer()
solution = scipy.integrate.solve_ivp(drdt_withks, t_span=(0,15000),
y0=(l10, l130, l150, l1350, c0, s0, cs0, l20, sa0, sb0, p0), method='RK45', rtol=1e-6)
end=timer()

#create list of all species from solved ODEs
t_ode45= solution.t
l1_ode45= solution.y[0]
l13_ode45= solution.y[1]
l15_ode45= solution.y[2]
l135_ode45= solution.y[3]
c_ode45= solution.y[4]
s_ode45= solution.y[5]
cs_ode45= solution.y[6]
l2_ode45= solution.y[7]
sa_ode45= solution.y[8]
sb_ode45= solution.y[9]
p_ode45= solution.y[10]


##########
###plot###
##########
#whole
fig = plt.figure(figsize=(24, 8))
ax = fig.add_subplot(1,1,1)
ax.plot(t_ode45, (l1_ode45), color = 'green',label='pre-circularised linear RNA number', linestyle='-', linewidth=2)
ax.plot(t_ode45, (l13_ode45), color = 'lime',label="3' processed linear RNA", linestyle=':', linewidth=2)
ax.plot(t_ode45, (l15_ode45), color = 'darkgreen',label="5' processed linear RNA", linestyle=':', linewidth=2)
ax.plot(t_ode45, (l135_ode45), color = 'dodgerblue',label="3' & 5' processed linear RNA", linestyle=':', linewidth=2)
ax.plot(t_ode45, (c_ode45), color = 'blue',label='circular RNA number', linestyle='-', linewidth=2)
ax.plot(t_ode45, (s_ode45), color = 'red',label="siRNA", linestyle='--', linewidth=2)
ax.plot(t_ode45, (cs_ode45), color = 'purple',label="siRNA:cRNA complex", linestyle='--', linewidth=2)
ax.plot(t_ode45, (l2_ode45), color = 'cyan',label="linearised RNA", linestyle='-', linewidth=2)
ax.plot(t_ode45, (sa_ode45), color = 'orangered',label="cleaved RNA product 'a' ", linestyle='--', linewidth=2)
ax.plot(t_ode45, (sb_ode45), color = 'orange',label="cleaved RNA product 'b' ", linestyle='--', linewidth=2)
ax.plot(t_ode45, (p_ode45), color = 'darkgoldenrod',label='protein number', linestyle='-', linewidth=2)
ax.set_title('Species number with time: whole system', fontsize=30)
ax.set_xlabel('time', fontsize=30)
ax.set_ylabel('copy number', fontsize=30)
plt.tick_params(axis='y', labelsize=0, length = 0)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=30)
plt.savefig('Deterministic.png', dpi=300)  
