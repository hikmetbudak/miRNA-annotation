#!/usr/bin/python

# Copyright(c) 2016 - Hikmet Budak <hikmet.budak@montana.edu>, HaliseBusra <hbusra@sabanciuniv.edu>

# SUmirLocate.py - a script that uses SUmirScreen.py outputs to count the number 
# of occurences of putative miRNAs in a given genome (or any other fasta file) 

# mirna representation (IDs are in forms of miR3157 etc.)
# follow the instructions on line 110 to change miRNA IDs to the forms of miR3157-5 etc. 
#################################################

import sys

#################################################
usage = "SUmirLocate.py genome.fasta SUmirScreen.output.out.tbl"
#################################################

try: script,genome,filename=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

#--------------------
def alllocations(alist, x):
    i=0
    while True:
        i = alist.find(x,i)
        if i==-1: return
        yield i
        i+=1
#----------------------

with open(genome,'r') as f:content=f.read()
lines = content.split("\n")

geneID = []
geneSeq = []
done = 0
seq = ''
for line in lines:
    if line.startswith('>'):
        if seq != '':
            geneSeq.append(seq)
            seq = ''
        geneID.append(line)
        done = 1
    else:
        if done == 1:
            seq += line
geneSeq.append(seq)

#--------------------

#filename = "plantmiRNAs.txt.fsa.results.tbl.hairpins.tbl.out.tbl"
with open(filename,'r') as f:content=f.read()
lines = content.split("\n")

out1 = open(filename+'.edited','w')
out3=open(filename+'.expression.tbl','w')
senseSeq = {}
antiSeq = {}
mirnaID = []
mirnaIndex = {}
howmany = {}

for line in lines:
    if line.startswith('miR'):
        mirna = line.split()[0]
        mirna2=line.split()[-1]
        if mirna in mirnaID:
            a = int(mirnaIndex[mirna])+1
            mirnaIndex[mirna] = str(a)
        else:
            mirnaIndex[mirna] = '1'
            mirnaID.append(mirna)
        
        index = line.split()[3]

        seq = line.split()[2]
        newseq = ''
        newseq2 = ''
        for s in seq:
            if s =='U':
                newseq += 'A'
                newseq2 += 'T'
            elif s == 'A':
                newseq += 'T'
                newseq2 += s
            elif s == 'G':
                newseq += 'C'
                newseq2 += s
            elif s == 'C':
                newseq += 'G'
                newseq2 += s
            elif s == 'N':
                newseq += 'N'
                newseq2 += s
            else:
                print 'something is wrong with ' + seq
        newseq = newseq[::-1]

        if newseq2 in senseSeq:senseSeq[newseq2] += ('\t'+mirna+'-'+mirnaIndex[mirna]+','+index)
        else:senseSeq[newseq2] = str(mirna+'-'+mirnaIndex[mirna]+','+index)

        if newseq in antiSeq:antiSeq[newseq] += ('\t'+mirna+'-'+mirnaIndex[mirna]+','+index)
        else:antiSeq[newseq] = str(mirna+'-'+mirnaIndex[mirna]+','+index)
        
        summary = '\t'.join(line.split()[1:])
        for mir2 in mirna2.split(','):out3.write(mirna+'-'+mirnaIndex[mirna]+'\t'+line.split()[1]+'\t'+mir2+'\n')
    
        out1.write(mirna+'-'+mirnaIndex[mirna]+'\t'+summary+'\n')

#   If you want to count mirna Isomers instead of general mirnaIDs 
#    uncomment the 1st line and comment the 2nd line below        
#        newmirna = mirna+'-'+mirnaIndex[mirna]
        newmirna = mirna
        if newmirna not in howmany:
            howmany[newmirna] = [0, 0]
out1.close()
out3.close()
#--------------------------

out2 = open('summary-premirna-locations.csv', 'w')
out2.write('mirnaID\tindex\treadID\tpremirna location\tpremirna length\tstrand')
try:
    for seq in senseSeq:
        geneCount = -1
        for gene in geneSeq:
            geneCount += 1
            n1=n2=0
            loc1 = list(alllocations(gene,seq))
            if loc1 != []:
                for loc in loc1:
                    allmirnas = senseSeq[seq].split('\t')
                    for mirnas in allmirnas:
                        mirna = '-'.join(mirnas.split(',')[0].split('-')[:-1])
                        index = mirnas.split(',')[1]
                        out2.write('\n'+mirnas.split(',')[0]+'\t'+index+'\t'+geneID[geneCount].split()[0][1:]+'\t'+str(loc)+'\t'+str(len(seq))+'\tSENSE')
                        howmany[mirna][0] += 1

    for seq in antiSeq:
        geneCount = -1
        for gene in geneSeq:
            geneCount += 1
            n1=n2=0
            loc1 = list(alllocations(gene,seq))
            if loc1 != []:
                for loc in loc1:
                    allmirnas = antiSeq[seq].split('\t')
                    for mirnas in allmirnas:
                        mirna = '-'.join(mirnas.split(',')[0].split('-')[:-1])
                        index = mirnas.split(',')[1]
                        out2.write('\n'+mirnas.split(',')[0]+'\t'+index+'\t'+geneID[geneCount].split()[0][1:]+'\t'+str(loc)+'\t'+str(len(seq))+'\tANTISENSE')
                        howmany[mirna][1] += 1
except:
    print "Check if there is only one miRNA per line.\nYou should choose one from the miRNAs separated by commas.\nError in mirnaID "+str(mirnas)
out2.close()

#---------------------------------------

output3 = open('summary-premirna-counts.csv','w')
output3.write('mirnaID\ton SENSE\ton ANTISENSE\ttotal')
for many in howmany:
    output3.write('\n%s\t%d\t%d\t%d' % (many,howmany[many][0],howmany[many][1],howmany[many][0]+howmany[many][1]))
output3.close()

