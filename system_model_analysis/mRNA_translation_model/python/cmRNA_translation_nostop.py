import numpy as np
import matplotlib.pyplot as plt 
import time
#set random seed
np.random.seed(1)

#function to specifu if number is mutliple of length of protein
def isMultipleof186(n): #aa in: gfp(140) + ribozyme linker regions(16) + 2a self cleaving peptide(20)
      
    while ( n > 0 ): 
        n = n - 186
  
    if ( n == 0 ): 
        return 1
  
    return 0

#how long to run simulation for, stays outside function as refer later to calulate mean
timespan = 40000
#gillespie function
def gillespie(kc, kd, ktra, kd4, c, rib, crib, aa, p, truncp, t):
    #output lists, each lsit contains species counts
    clist = [1]
    riblist = [1]
    criblist = [0]
    aalist = [0]
    plist = [0]
    truncplist = [0]
    truncplenlist = []
    tlist = [0] 
    #output is a list of each species list
    output = [tlist,clist,riblist,criblist,plist,truncplist,truncplenlist]
    #get species count until reach end of timespan

    while t < timespan:  #                
        #make empty list of reaction propensities
        rlist = []
        #define reactions and what they do
        reactiona = kc*c*rib # crib+1, c-1, rib-1 ribosome associaiton
        reactionb = kd*crib # crib-1, c+1, rib+1 ribosome dissocaition
        reactionc = ktra*crib #  aa+1 transaltion elongation
        reactiond = kd4*p # p-1 protein degradation
        reactione = kd4*truncp # p-1 truncated prtoein degradation
        #extend empty rlist with reaction propensities
        rlist.extend(value for name, value in sorted(locals().items(), key=lambda item: item[0]) if name.startswith('reaction'))
        #calculate sum of propensities
        R = sum(rlist)
        #propensity divided by sum of propensities, so get probabiltity of reaction occuring
        xlist = [x / R for x in rlist]
        #get ranges of probabilties, so if random no. falls in a range that reaction happens
        problist = [sum(xlist[0:i+1]) for i,v in enumerate(xlist)]
        #what letter appends to what list
        clist.append(c)
        riblist.append(rib)
        criblist.append(crib)
        aalist.append(aa)
        plist.append(p)
        truncplist.append(truncp)
        #pick random no. between 0,1
        rand = np.random.rand()
        #uncomment these if want real time update on species
        # print("aalist: ", aalist)
        # print("clist: ", clist)
        # print("riblist: ", riblist)
        # print("criblist: ", criblist)
        # print("plist: ", p)
        # print("truncplist: ", plist)
        # print("c: ", c)
        # print("rib: ", rib)
        # print("crib: ", crib)
        # print("aa: ", aa)
        # print("p: ", p)
        # print("truncp: ", truncp)
        # print("aa: ", aa)
        # print("kc: ", kc)
        # print("kd: ", kd)
        # print("ktra: ", ktra)
        # print("kd4: ", kd4)

        # print("length of protein: ", len(aalist))
        #if falls between x1 and x2, add/minus 1 from list/lists
        #if aa count is mutliple of length of protein add +1 to equivalent protein list
        if isMultipleof186(aa):   
            p = plist[-1] + 1
            crib= criblist[-1] + 0 
            c = clist[-1] + 0                     
            rib= riblist[-1] + 0                                     
        if rand >=0 and rand <= problist[0]: #reactiona 
            c = clist[-1] - 1                     
            rib= riblist[-1] - 1                   
            crib= criblist[-1] + 1   
            print("reactiona: ribosome:rna association")
        elif rand >problist[0] and rand <= problist[1]: #reactionb
            c = clist[-1] + 1                     
            rib= riblist[-1] + 1                   
            crib= criblist[-1] - 1  
            truncp = truncplist[-1]+1 #add number of polyproteins to lsit
            truncplenlist.append(len(aalist)) #add length of protein to list
            aalist.clear() 
            aalist = [0]  
            print("reactionb: spontaneous ribosome dissociation")
        elif rand >problist[1] and rand <= problist[2]: #reactionc
            aa = aalist[-1] + 1
            print("reactionc: translation elongation")
        elif rand >problist[2] and rand <= problist[3]: #reactione  
            p = plist[-1] - 1
            c = clist[-1] + 0
            rib = riblist[-1] + 0
            crib = criblist[-1] + 0
            print("reactiond: protein degradation")
        elif rand >problist[3] and rand <= problist[4]: #reactione  
            truncp = truncplist[-1] - 1   
            c = clist[-1] + 0
            rib = riblist[-1] + 0
            crib = criblist[-1] + 0
            print("reactione: truncated protein degradation")        
        #time increments
        tau = np.random.exponential(1/R)
        t = tlist[-1] + tau 
        tlist.append(t)
        print("time elapsed:", t)
        print("___________________________________")
    return output                              



kc = 0.333 # tir estiamted at 20 per minute
kd = 0.003 # 2.4x10-4 per codon, 140 codons, will take 9.333s to translate whole protein |0.028/9.333|
ktra = 15 # 15 aa per second elongation
kd4 = 1/(110*60) #1/(110*60) # unstable gfp |110 minute hl|
print(kd4)
#initial species counts
t = 0
c = 1
rib = 1
crib = 0
aa = 0
p = 0    
truncp = 0

#no. of simualtions
no_of_sim = 5
simulations = []

#run simualtions
time_start = time.time()
for x in range(no_of_sim):
    simulations = simulations + gillespie(kc, kd, ktra, kd4, c, rib, crib, aa, p, truncp, t)
time_end = time.time()
print("time elapse: ", time_start - time_end)

f = int(len(simulations)/no_of_sim)


#calcualting steady state of each species
t_range = np.arange(0, timespan)

#get arrays of lsits of time, equivalent proteins and polyprotiens
t_array = [np.array(x) for x in simulations[::f]]
p_array = [np.array(x) for x in simulations[4::f]]
tp_array = [np.array(x) for x in simulations[5::f]]

#get average of copy numbers equivalents and polyproeins
p_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(p_array,t_array)]) for t in t_range ]
tp_average_list = [np.mean([x[(d<=t).sum()] for x,d in zip(tp_array,t_array)]) for t in t_range ]

print("the steady state equivalent of: ", p_average_list[-1], " are translated")
print("steady state of polyproteins is: ", tp_average_list[-1])

###########
###plots###
###########
counter = 0
plt.figure(figsize=(12, 8))
while counter <= ((no_of_sim*f)-1):
#    plt.plot(simulations[counter],simulations[counter+1], label = "circular RNA" if counter == 0 else "", lw = 1, color = 'blue', linestyle=":") #l1
#    plt.plot(simulations[counter],simulations[counter+2], label = "ribosome" if counter == 0 else "", lw = 1, color = "red", linestyle=":") #c
#    plt.plot(simulations[counter],simulations[counter+3], label = "circular RNA:ribosome" if counter == 0 else "", lw = 1, color = 'purple') #l1
    plt.plot(simulations[counter],simulations[counter+4], label = "equivalent single unit protein" if counter == 0 else "", lw = 1, color = "green") #c    
    plt.plot(simulations[counter],simulations[counter+5], label = "polyprotein" if counter == 0 else "", lw = 1, color = "green", linestyle=":") #c    
    plt.xlabel('Time (s)', fontsize=20)
    plt.ylabel('copy no.', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.grid()
    plt.legend(fontsize = 20)
    counter += f
plt.show()

#get hisotrgam of size of polyproteins
plt.hist(simulations[+6], density=False, bins=100, color="dodgerblue")  # `density=False` would make counts
plt.hist(simulations[+13], density=False, bins=100, color="dodgerblue")  # `density=False` would make counts
plt.hist(simulations[+20], density=False, bins=100, color="dodgerblue")  # `density=False` would make counts
plt.hist(simulations[+27], density=False, bins=100, color="dodgerblue")  # `density=False` would make counts
plt.hist(simulations[+34], density=False, bins=100, color="dodgerblue")  # `density=False` would make counts
plt.xticks(np.arange(0, 45000, step=1500), rotation = "vertical")
plt.ylabel('Frequency')
plt.xlabel('Polyrotein length (aa)')
plt.show()


###

#TO GET RIBOSOME:MRNA ASSOCIATION GRAPHS, SET timespan TO 1000 and uncomment code below and comment out graph code above
#AND SET no_of_sim TO 1
# counter = 0
# plt.figure(figsize=(12, 8))
# while counter <= ((no_of_sim*f)-1):
#     plt.plot(simulations[counter],simulations[counter+1], label = "circular RNA" if counter == 0 else "", lw = 2, color = 'blue', linestyle="-") #l1
#     plt.plot(simulations[counter],simulations[counter+2], label = "ribosome" if counter == 0 else "", lw = 2, color = "red", linestyle=":") #c
# #    plt.plot(simulations[counter],simulations[counter+3], label = "circular RNA:ribosome" if counter == 0 else "", lw = 2, color = 'purple') #l1
# #     plt.plot(simulations[counter],simulations[counter+4], label = "protein" if counter == 0 else "", lw = 1, color = "green") #c    
# #     plt.plot(simulations[counter],simulations[counter+5], label = "truncated protein" if counter == 0 else "", lw = 1, color = "green", linestyle=":") #c    
#     plt.xlabel('Time (s)', fontsize=20)
#     plt.ylabel('copy no.', fontsize=20)
#     plt.grid()
#     plt.legend(fontsize = 20)
#     plt.tick_params(axis='both', which='major', labelsize=15)
#     counter += f
# plt.show()


# counter = 0
# plt.figure(figsize=(25, 8))
# while counter <= ((no_of_sim*f)-1):
# #    plt.plot(simulations[counter],simulations[counter+1], label = "circular RNA" if counter == 0 else "", lw = 2, color = 'blue', linestyle=":") #l1
# #    plt.plot(simulations[counter],simulations[counter+2], label = "ribosome" if counter == 0 else "", lw = 2, color = "red", linestyle=":") #c
#     plt.plot(simulations[counter],simulations[counter+3], label = "circular RNA:ribosome" if counter == 0 else "", lw = 2, color = 'purple') #l1
# #     plt.plot(simulations[counter],simulations[counter+4], label = "protein" if counter == 0 else "", lw = 1, color = "green") #c    
# #     plt.plot(simulations[counter],simulations[counter+5], label = "truncated protein" if counter == 0 else "", lw = 1, color = "green", linestyle=":") #c    
#     plt.xlabel('Time (s)', fontsize=20)
#     plt.ylabel('copy no.', fontsize=20)
#     plt.grid()
#     plt.legend()
#     counter += f
# plt.show()
