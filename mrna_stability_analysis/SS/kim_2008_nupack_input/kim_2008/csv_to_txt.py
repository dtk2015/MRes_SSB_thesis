#this code creates individual .IN files for nput into NUAPCK
#the code also creates shell scripts that can be use NUPACK on all .IN fiels in linux server

import csv

#to create files to get mfe structure
mfe = []

with open('kim_2008_seq.csv') as csvfile:
    
    reader = csv.DictReader(csvfile)
    for row in reader:
        file_name ='{0}.in'.format(row['Gene ID']) 
        mfe.append(file_name)
        print(len(mfe))
        with open(file_name, 'w') as f:
            f.write(row["seq_UTR+25cds"]) 

            
mfe = [s.strip('.in') for s in mfe]

            


runjobspfunc= open("runjobspfunc","w+")
runjobspfunc.write("\n")
runjobspfunc.write("#############################################")
runjobspfunc.write("\n")
runjobspfunc.write("##### get partition function ####################")
runjobspfunc.write("\n")
runjobspfunc.write("#############################################")
runjobspfunc.write("\n")

for i in mfe:
    p = "pfunc -T 37 " + "/home/dp2015/kim_2008/mfe/"
    p += i
    p += " "    
    p += "> " + "/home/dp2015/kim_2008/output/"
    p += i + ".pfunc"
    runjobspfunc.write(p)
    runjobspfunc.write(" ;")
    runjobspfunc.write("\n")

runjobsmfe= open("runjobsmfe","w+")
runjobsmfe.write("#############################################")
runjobsmfe.write("\n")
runjobsmfe.write("##### get mfe structures ####################")
runjobsmfe.write("\n")
runjobsmfe.write("#############################################")
runjobsmfe.write("\n")

for i in mfe:
    p = "mfe -T 37 " + " "+"/home/dp2015/kim_2008/mfe/"
    p += i    
    runjobsmfe.write(p)
    runjobsmfe.write(" ;")
    runjobsmfe.write("\n")

runjobsppairs= open("runjobsppairs","w+")
runjobsppairs.write("\n")
runjobsppairs.write("#############################################")
runjobsppairs.write("\n")
runjobsppairs.write("##### get pairs #############################")
runjobsppairs.write("\n")
runjobsppairs.write("#############################################")
runjobsppairs.write("\n")

for i in mfe:
    p = "pairs -T 37 " + " "+"/home/dp2015/kim_2008/mfe/"
    p += i    
    runjobsppairs.write(p)
    runjobsppairs.write(" ;")
    runjobsppairs.write("\n")



for i in mfe:
    p = "mv" + " "+"/home/dp2015/kim_2008/mfe/"
    p += i +".mfe"   
    p += " " + "/home/dp2015/kim_2008/output/" + i + ".mfe"
    runjobsmfe.write(p)
    runjobsmfe.write(" ;")
    runjobsmfe.write("\n")

# runjobsppairs.write("\n")
# runjobsppairs.write("#############################################")
# runjobsppairs.write("\n")
# runjobsppairs.write("##### move results to output folder ########")
# runjobsppairs.write("\n")
# runjobsppairs.write("#############################################")
# runjobsppairs.write("\n")

# for i in mfe5:
#     p = "mv" + " "+"/home/dp2015/kim_2008/mfe/"
#     p += i +".ppairs"   
#     p += " " + "/home/dp2015/kim_2008/output/" + i + ".ppairs"
#     runjobsppairs.write(p)
#     runjobsppairs.write(" ;")
#     runjobsppairs.write("\n")

# runjobssubopt.write("\n")
# runjobssubopt.write("#############################################")
# runjobssubopt.write("\n")
# runjobssubopt.write("##### move results to output folder ########")
# runjobssubopt.write("\n")
# runjobssubopt.write("#############################################")
# runjobssubopt.write("\n")

# for i in sub5:
#     p = "mv" + " "+"/home/dp2015/kim_2008/sub/"
#     p += i +".subopt"   
#     p += " " + "/home/dp2015/kim_2008/output/" + i + ".subopt"
#     runjobssubopt.write(p)
#     runjobssubopt.write(" ;")
#     runjobssubopt.write("\n")









