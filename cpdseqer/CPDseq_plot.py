#!/usr/bin/python
import gzip
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
from datetime import datetime

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-d", "--dinucleotide", help = "input dinucleotide count file list(Tab delimited)", required=True)
parser.add_argument("-b", "--interval", help = "input interval bed file", required=True)
#parser.add_argument("-o", "--output", help = "output file for the results (required)", required=True)
args = parser.parse_args()


def mutationCount(dinucleotide,interval,max_interval):
    chr_last="0"
    mutation=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    cnt=pd.DataFrame(0, index=range(1,max_interval+1), columns=mutation)
    D={}
    for line1 in interval:
        tokens1=line1.split("\t")
        chr=tokens1[0]
        if chr != chr_last: 
            if len(D) != 0:
                time = datetime.now()
                print time.strftime("%D %H:%M:%S") + " - Processing chromosome "+ chr_last +"..."
                dinucleotide.seek(0)
                for line2 in dinucleotide:
                    tokens2=line2.split("\t")
                    key=tokens2[0]+"_"+tokens2[1]
                    if key in D and tokens2[3] in mutation:
                        cnt.ix[D[key],tokens2[3]] += 1
                D.clear()
            chr_last=chr
        start=int(tokens1[1])
        end=int(tokens1[2])
        k=1
        for pos in range(start, end+1):
            key=tokens1[0]+"_"+str(pos)
            D[key]=k
            k=k+1
    if len(D) != 0:
        time = datetime.now()
        print time.strftime("%D %H:%M:%S") + " - Processing chromosome "+ chr_last +"..."
        dinucleotide.seek(0)
        for line2 in dinucleotide:
            tokens2=line2.split("\t")
            key=tokens2[0]+"_"+tokens2[1]
            if key in D:
                cnt.ix[D[key],tokens2[3]] += 1
    D.clear()
    return cnt
    

##interval
print "Accessing chromosome information of " + args.interval + " file and performing sort..."
GZIP_MAGIC_NUMBER = "1f8b"
f = open(args.interval)
if f.read(2).encode("hex")==GZIP_MAGIC_NUMBER:
        command="zcat "+args.interval+" | sort -V -k1,1 -k2,2 > temp_interval"
else:
        command="sort -V -k1,1 -k2,2 "+args.interval+" >temp_interval"
os.system(command)
f.close()

interval=open("temp_interval", "r")
last_chr="0"
last_end=0
overlap_counter=0
max_interval=0
n=0
for line in interval:
    tokens=line.split("\t")
    n += 1
    if len(tokens)<3:
        sys.exit("\n...\n Format error: " + line + " \n The interval file should contain at least three columns, each representing chromosome,start position and end position, such as: \n chr1	826186	826332 \n...\n")
    if tokens[0]==last_chr and int(tokens[1])<=last_end:
        overlap_counter += 1
    last_chr=tokens[0]
    last_end=int(tokens[2])
    cur_interval=int(tokens[2])-int(tokens[1])+1
    if cur_interval>max_interval:
        max_interval=cur_interval
interval.seek(0)
if overlap_counter>0:
    print "\n...\n WARNING: there are " + str(overlap_counter) + " out of " + str(n) + " overlapping intervals in "+ args.interval +"\n...\n"

##dinucleotide
list_dinucleotide=open(args.dinucleotide)
for line in list_dinucleotide:
    tokens=line.split("\t")
    f=open(tokens[0]) 
    if f.read(2).encode("hex")==GZIP_MAGIC_NUMBER:
            dinucleotide=gzip.open(tokens[0], "r")
    else:
            dinucleotide=open(tokens[0], "r")
    f.close()
    print "\n Processing dinucleotide file: "+ tokens[0]
    interval.seek(0)
    cnt=mutationCount(dinucleotide,interval,max_interval)
    cnt.to_csv(os.path.splitext(tokens[0])[0]+".mutation.count.txt", sep='	', index=True, header=True)
    dinucleotide.close()
interval.close()

command="rm temp_interval"
os.system(command)
                                
