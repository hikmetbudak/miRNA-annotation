#!/usr/bin/python

# Copyright(c) 2016 - Hikmet Budak <hikmet.budak@montana.edu>, HaliseBusra <hbusra@sabanciuniv.edu>

# SUmirScreen.py - a script that uses SUmirFold.pl outputs to screen structures of predicted miRNAs for specified folding criteria

# folding criteria:
# 1- Mismatches, max 4 for sense and max 6 for antisense strands
# 2- Dicer cut points (mature mirna start and mirna* start sites)
# 3- Multiloop (only one loop is acceptable between mature mirna and mirna*)
# 4- Head (both mature mirna and mirna* sequences should not be involved in the head part)
#################################################

import sys
import os

#################################################
usage = "SUmirScreen.py SUmirFold.hairpins.tbl.fileoutput SUmirFold.hairpins.folderoutput"
#################################################

try: script,filename,dic=sys.argv
except:
    try: script,filename=sys.argv
    except: sys.exit("Correct usage is:\n"+usage)
    dic=filename[:-4]
    if not os.path.isdir(dic): sys.exit("Correct usage is:\n"+usage)

file = open(filename,'r')
content = file.read()
file.close()

lines = content.split("\n")
lineID = -1
hairpins = []
start1 = []
end1 = []
screen = []
whole_data = []
blacklist = ['Unique','Hit']
try:
    x = -1
    for line in lines:
        lineID += 1

        if not line.startswith('\t') and line != '' and line.split()[0] not in blacklist:
            x += 1
            store = line.split("\t")
            n = len(store)
            temp = line
            
            if n<20:
                temp += lines[lineID+1]
                n += len(lines[lineID+1].split('\t'))
                if n < 20:
                    temp += lines[lineID+2]
            store = temp.split('\t')
            whole_data.append(temp)
            hairpins.append(dic + "/" + store[0] + ".hairpin.fsa_1.ct")
            start1.append(int(store[8]))
            end1.append(int(store[9]))
            screen.append('')
except:
    print '1: error in file: ' + filename + '\t in line: ' + str(line)

er2count=0
count = -1
for i in hairpins:
    count += 1
    
    file = open(i,'r')
    content = file.read()
    file.close()
    
    lines = content.split('\n')
    
# store a and b strands
    a = [0]
    b = [0]
    c = ['0']
    m = 0
    try:
        for line in lines:
            m+= 1
            if m>1 and line != '':
                store = line.split('\t')
                a.append(int(store[0]))
                b.append(int(store[4]))
                c.append(store[1])
            elif m ==1:
                primirnalength = int(line.split('\t')[0])
    
        cont = True    
        while cont:
                    
# 1- check for mismatches, max 4 for sense and max 6 for antisense
            mismatchB = 0
            start2 = b[end1[count]-2]
            end2 = b[start1[count]]+2
            mirna2=''
            for j in range(start2,end2+1):
                mirna2+=c[j]
                if b[j] == 0:
                    mismatchB += 1
                    
                    if mismatchB > 6:
                        screen[count] = 'MismatchB'
                        break
            newdata=whole_data[count].split('\t')
            newdata[10]=str(start2)
            newdata[11]=str(end2)
            newdata[12]=mirna2
            whole_data[count]='\t'.join(newdata)
            if screen[count] == 'MismatchB':
                break
            
            mismatchA = 0
            for j in range(start1[count],end1[count]+1):
                if b[j] == 0:
                    mismatchA += 1
                
                    if mismatchA > 4:
                        screen[count] = 'Mismatch'
                        break
    
            if screen[count] == 'Mismatch': break
       
# 2- check for Dicer functions
            if b[start1[count]] == 0:
                screen[count] = 'Dicer'
                break
            if b[end1[count]-2] == 0:
                screen[count] = 'Dicer'
                break
    
# determine start-end sites
            start = start1[count]
            end = end2
            if start2 < start:
                start = start2
                end = end1[count]

# 3- check for multiloop
            temp = b[start]
            for k in range(start+1,end+1):
                check = int(b[k])
                if check != 0:
                    if check < temp:
                        temp = check
                    elif check > temp:
                        screen[count] = 'Multiloop'
                        break

            if screen[count] == 'Multiloop': break

# 4- check for head                
            temp = b[start]
            check = 0
            for k in range(start+1,end+1):
                if b[k] != 0:
                    if check == b[k]:
                        u2 = []
                        for u in range(k-x,k): u2.append(u)
                        for test in range(start1[count],end1[count]+1):
                            if test in u2:
                                screen[count] = 'Head'
                                break
                        for test in range(start2,end2+1):
                            if test in u2:
                                screen[count] = 'Head'
                                break
                    x = 0
                elif b[k] == 0:
                    x += 1
                    if x == 1: check = a[k-1]

            if screen[count] == 'Head': break
                    
            screen[count] = 'OK'
            break
    except:
        er2count+=1
        #print '2: error in file: ' + i
        
output = open((filename + '.edited.tbl'),'w')

starters = 'Screen\tUnique Hit ID\tNew miRNA ID\tNew miRNA Sequence\tNew miRNA Length\tConserved miRNA ID\tConserved miRNA Sequence\tConserved miRNA Mismatch\tSequence ID\tMature Start\tMature End\tmiRNA* Start\tmiRNA* End\tmiRNA* Sequence\tHairpin Location\tPre-miRNA Length\tPre-miRNA MFE\tPre-miRNA GC%\tPre-miRNA MFEI\tPre-miRNA Start\tPre-miRNA Sequence\tmirna* Sequence'
output.write(starters)

mirnas = []
seqs = []
mismatch = []
homolog=[]
index = []
location = []
mirna2=[]
for i in range(0,len(screen)):
    data = whole_data[i]
    output.write('\n' + screen[i] +'\t' + data)
    start2 = b[end1[count]-2]
    end2 = b[start1[count]]+2
    if screen[i] == 'OK' and not ('N' in data.split('\t')[2]):
        try:
            seqs.append(data.split('\t')[2] + '\t' + data.split('\t')[19])
            mirname=data.split('\t')[1][10:]
            if mirname.startswith("miR-"):
                mirID = mirname.split('-')[1].split('.')[0]
                mir = "miR-" + ''.join(c for c in mirID if c in '0123456789')
            elif mirname.startswith("bantam"):
                mir="bantam"
            else: 
                mirID = mirname.split('-')[0].split('.')[0]
                mir = "miR"+''.join(c for c in mirID if c in '0123456789')
            mismatch.append(data.split('\t')[6])
            homolog.append(data.split('\t')[4])
            index.append(str(i+1))
            loc=data.split('\t')[13][:1]
            mirna2.append(data.split('\t')[12])
            location.append(loc)
            homolog_loc = data.split('\t')[4].split('-')[-1]
            if (homolog_loc in ['3p','5p']) and loc!=homolog_loc[:1]: mir = mir
            else: mir += '-%sp'%(loc)
            mirnas.append(mir)
            
        except:
            print '3: error in: ' + data
output.close()

full = {}
for i in range(0,len(mirnas)):
    if seqs[i] in full:
        if full[seqs[i]][1] > mismatch[i]:
            full[seqs[i]] = [mirnas[i], mismatch[i], index[i],homolog[i],mirna2[i]]
        elif full[seqs[i]][1] == mismatch[i]:
            if mirnas[i] not in full[seqs[i]][0].split(','):
                full[seqs[i]][0] += ',' + mirnas[i]
                full[seqs[i]][2] += ',' + index[i]
                full[seqs[i]][3] += ',' + homolog[i]
                if mirna2[i] not in full[seqs[i]][4].split(','):
                    full[seqs[i]][3] += ',' + mirna2[i]
                print 'similar mismatch found. look at mirnaID: ' + full[seqs[i]][0] 
    else:
        full[seqs[i]] = [mirnas[i], mismatch[i], index[i],homolog[i],mirna2[i]]

output2 = open((filename + '.out.tbl'),'w')
for x in full: output2.write(full[x][0] + '\t' + x + '\t' + '\t'.join(full[x][2:])+'\n')
output2.close()

print "Number of structures with missing information: "+str(er2count)
print "\nFinished successfully."