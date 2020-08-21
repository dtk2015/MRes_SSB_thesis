#this code creates individual .IN files for nput into NUAPCK
#the code also creates shell scripts that can be use NUPACK on all .IN fiels in linux server

import csv

#to create files to get mfe structure
mfe5 = []

with open('5utr_seq.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.in'.format("5" + row['mRNA_name']) 
        mfe5.append(file_name)
        with open(file_name, 'w') as f:
            f.write(row["5' UTR Sequence"])  
            
#to create fiels to get sub optiaml sturctures
sub5 = []

with open('5utr_seq.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.in'.format("subopt" + "5" + row['mRNA_name']) 
        sub5.append(file_name)
        with open(file_name, 'w') as f:
            f.write(row["5' UTR Sequence"])        
            f.write("\n")
            f.write(row["subopt"])

mfe5 = [s.strip('.in') for s in mfe5]
sub5 = [s.strip('.in') for s in sub5]

runjobs5pfunc= open("runjobs5pfunc","w+")
runjobs5pfunc.write("\n")
runjobs5pfunc.write("#############################################")
runjobs5pfunc.write("\n")
runjobs5pfunc.write("##### get partition function ####################")
runjobs5pfunc.write("\n")
runjobs5pfunc.write("#############################################")
runjobs5pfunc.write("\n")

for i in mfe5:
    p = "pfunc -T 37 " + "/home/dp2015/mRNA/5/mfe_5/"
    p += i
    p += " "    
    p += "> " + "/home/dp2015/mRNA/5/output/"
    p += i + ".pfunc"
    runjobs5pfunc.write(p)
    runjobs5pfunc.write(" ;")
    runjobs5pfunc.write("\n")

runjobs5mfe= open("runjobs5mfe","w+")
runjobs5mfe.write("#############################################")
runjobs5mfe.write("\n")
runjobs5mfe.write("##### get mfe structures ####################")
runjobs5mfe.write("\n")
runjobs5mfe.write("#############################################")
runjobs5mfe.write("\n")

for i in mfe5:
    p = "mfe -T 37" + " "+"/home/dp2015/mRNA/5/mfe_5/"
    p += i    
    runjobs5mfe.write(p)
    runjobs5mfe.write(" ;")
    runjobs5mfe.write("\n")

runjobs5ppairs= open("runjobs5ppairs","w+")
runjobs5ppairs.write("\n")
runjobs5ppairs.write("#############################################")
runjobs5ppairs.write("\n")
runjobs5ppairs.write("##### get pairs #############################")
runjobs5ppairs.write("\n")
runjobs5ppairs.write("#############################################")
runjobs5ppairs.write("\n")

for i in mfe5:
    p = "pairs -T 37" + " "+"/home/dp2015/mRNA/5/mfe_5/"
    p += i    
    runjobs5ppairs.write(p)
    runjobs5ppairs.write(" ;")
    runjobs5ppairs.write("\n")

runjobs5subopt= open("runjobs5subopt","w+")
runjobs5subopt.write("\n")
runjobs5subopt.write("#############################################")
runjobs5subopt.write("\n")
runjobs5subopt.write("##### get sub optimal mfe structures ########")
runjobs5subopt.write("\n")
runjobs5subopt.write("#############################################")
runjobs5subopt.write("\n")


for i in sub5:
    p = "subopt -T 37e" + " "+"/home/dp2015/mRNA/5/sub_5/"
    p += i    
    runjobs5subopt.write(p)
    runjobs5subopt.write(" ;")
    runjobs5subopt.write("\n")

runjobs5mfe.write("\n")
runjobs5mfe.write("#############################################")
runjobs5mfe.write("\n")
runjobs5mfe.write("##### move results to output folder ########")
runjobs5mfe.write("\n")
runjobs5mfe.write("#############################################")
runjobs5mfe.write("\n")

for i in mfe5:
    p = "mv" + " "+"/home/dp2015/mRNA/5/mfe_5/"
    p += i +".mfe"   
    p += " " + "/home/dp2015/mRNA/5/output/" + i + ".mfe"
    runjobs5mfe.write(p)
    runjobs5mfe.write(" ;")
    runjobs5mfe.write("\n")

runjobs5ppairs.write("\n")
runjobs5ppairs.write("#############################################")
runjobs5ppairs.write("\n")
runjobs5ppairs.write("##### move results to output folder ########")
runjobs5ppairs.write("\n")
runjobs5ppairs.write("#############################################")
runjobs5ppairs.write("\n")

for i in mfe5:
    p = "mv" + " "+"/home/dp2015/mRNA/5/mfe_5/"
    p += i +".ppairs"   
    p += " " + "/home/dp2015/mRNA/5/output/" + i + ".ppairs"
    runjobs5ppairs.write(p)
    runjobs5ppairs.write(" ;")
    runjobs5ppairs.write("\n")

runjobs5subopt.write("\n")
runjobs5subopt.write("#############################################")
runjobs5subopt.write("\n")
runjobs5subopt.write("##### move results to output folder ########")
runjobs5subopt.write("\n")
runjobs5subopt.write("#############################################")
runjobs5subopt.write("\n")

for i in sub5:
    p = "mv" + " "+"/home/dp2015/mRNA/5/sub_5/"
    p += i +".subopt"   
    p += " " + "/home/dp2015/mRNA/5/output/" + i + ".subopt"
    runjobs5subopt.write(p)
    runjobs5subopt.write(" ;")
    runjobs5subopt.write("\n")
