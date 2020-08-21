#this code creates individual .IN files for nput into NUAPCK
#the code also creates shell scripts that can be use NUPACK on all .IN fiels in linux server

import csv

#to create files to get mfe structure
mfe3 = []
with open('3utr_seq.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.in'.format("3" + row['mRNA_name']) 
        mfe3.append(file_name)
        with open(file_name, 'w') as f:
            f.write(row["3' UTR Sequence"])  
#to create fiels to get sub optiaml sturctures
sub3 = []
with open('3utr_seq.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.in'.format("subopt" + "3" + row['mRNA_name']) 
        sub3.append(file_name)
        with open(file_name, 'w') as f:
            f.write(row["3' UTR Sequence"])        
            f.write("\n")
            f.write(row["subopt"])
            
mfe3 = [s.strip('.in') for s in mfe3]
sub3 = [s.strip('.in') for s in sub3]

runjobs3pfunc= open("runjobs3pfunc","w+")
runjobs3pfunc.write("\n")
runjobs3pfunc.write("#############################################")
runjobs3pfunc.write("\n")
runjobs3pfunc.write("##### get partition function ####################")
runjobs3pfunc.write("\n")
runjobs3pfunc.write("#############################################")
runjobs3pfunc.write("\n")

for i in mfe3:
    p = "pfunc -T 37 " + "/home/dp2015/mRNA/3/mfe_3/"
    p += i
    p += " "    
    p += "> " + "/home/dp2015/mRNA/3/output/"
    p += i + ".pfunc"
    runjobs3pfunc.write(p)
    runjobs3pfunc.write(" ;")
    runjobs3pfunc.write("\n")

runjobs3mfe= open("runjobs3mfe","w+")
runjobs3mfe.write("#############################################")
runjobs3mfe.write("\n")
runjobs3mfe.write("##### get mfe structures ####################")
runjobs3mfe.write("\n")
runjobs3mfe.write("#############################################")
runjobs3mfe.write("\n")

for i in mfe3:
    p = "mfe -T 37" + " "+"/home/dp2015/mRNA/3/mfe_3/"
    p += i    
    runjobs3mfe.write(p)
    runjobs3mfe.write(" ;")
    runjobs3mfe.write("\n")

runjobs3ppairs= open("runjobs3ppairs","w+")
runjobs3ppairs.write("\n")
runjobs3ppairs.write("#############################################")
runjobs3ppairs.write("\n")
runjobs3ppairs.write("##### get pairs #############################")
runjobs3ppairs.write("\n")
runjobs3ppairs.write("#############################################")
runjobs3ppairs.write("\n")
for i in mfe3:
    p = "pairs -T 37" + " "+"/home/dp2015/mRNA/3/mfe_3/"
    p += i    
    runjobs3ppairs.write(p)
    runjobs3ppairs.write(" ;")
    runjobs3ppairs.write("\n")


runjobs3subopt= open("runjobs3subopt","w+")
runjobs3subopt.write("\n")
runjobs3subopt.write("#############################################")
runjobs3subopt.write("\n")
runjobs3subopt.write("##### get sub optimal mfe structures ########")
runjobs3subopt.write("\n")
runjobs3subopt.write("#############################################")
runjobs3subopt.write("\n")


for i in sub3:
    p = "subopt -T 37" + " "+"/home/dp2015/mRNA/3/sub_3/"
    p += i    
    runjobs3subopt.write(p)
    runjobs3subopt.write(" ;")
    runjobs3subopt.write("\n")

runjobs3mfe.write("\n")
runjobs3mfe.write("#############################################")
runjobs3mfe.write("\n")
runjobs3mfe.write("##### move results to output folder ########")
runjobs3mfe.write("\n")
runjobs3mfe.write("#############################################")
runjobs3mfe.write("\n")

for i in mfe3:
    p = "mv" + " "+"/home/dp2015/mRNA/3/mfe_3/"
    p += i +".mfe"   
    p += " " + "/home/dp2015/mRNA/3/output/" + i + ".mfe"
    runjobs3mfe.write(p)
    runjobs3mfe.write(" ;")
    runjobs3mfe.write("\n")

runjobs3ppairs.write("\n")
runjobs3ppairs.write("#############################################")
runjobs3ppairs.write("\n")
runjobs3ppairs.write("##### move results to output folder ########")
runjobs3ppairs.write("\n")
runjobs3ppairs.write("#############################################")
runjobs3ppairs.write("\n")

for i in mfe3:
    p = "mv" + " "+"/home/dp2015/mRNA/3/mfe_3/"
    p += i +".ppairs"   
    p += " " + "/home/dp2015/mRNA/3/output/" + i + ".ppairs"
    runjobs3ppairs.write(p)
    runjobs3ppairs.write(" ;")
    runjobs3ppairs.write("\n")

runjobs3subopt.write("\n")
runjobs3subopt.write("#############################################")
runjobs3subopt.write("\n")
runjobs3subopt.write("##### move results to output folder ########")
runjobs3subopt.write("\n")
runjobs3subopt.write("#############################################")
runjobs3subopt.write("\n")


for i in sub3:
    p = "mv" + " "+"/home/dp2015/mRNA/3/sub_3/"
    p += i +".subopt"   
    p += " " + "/home/dp2015/mRNA/3/output/" + i + ".subopt"
    runjobs3subopt.write(p)
    runjobs3subopt.write(" ;")
    runjobs3subopt.write("\n")
