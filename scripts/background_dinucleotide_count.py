# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:57:32 2020

@author: lijiang
"""

from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-fa", "--inputfastaFile", help = "input fasta file name", required=True)
parser.add_argument("-bed", "--inputbedFile", help = "input bed file name", default='NONE')
parser.add_argument("-o", "--outputbedFile", help = "output file name", required=False)
args = parser.parse_args()


fasta_file = args.inputfastaFile

#fasta_file = "seq_unique.fasta"
seqs = list(SeqIO.parse(fasta_file, "fasta"))

lags = []
for ii in range(len(seqs)):
    lags.append(seqs[ii].description)


    
bedfile = args.inputbedFile
#bedfile = "test.bed"


result_all_all = pd.DataFrame(index = ['AA', 'AC', 'AT', 'AG', 'CC', 'CA', 'CT', 'CG', 'GG', 'GA', 'GC', 'GT','TT', 'TA', 'TC', 'TG'])


if bedfile != 'NONE':
    data_chrom = []
    start_pos = []
    end_pos = []
    types = []
    type_lag = "No"
    
    f = open(bedfile, 'r')    
    words = f.readline().split() 
    
    if len(words) == 4:        
        type_lag = "Yes"    
    
    if words[1].isdigit():
        words[-1] = words[-1].strip()             
        data_chrom.append(words[0])
        start_pos.append(int(words[1]))
        end_pos.append(int(words[2]))
        if len(words) == 4:        
            types.append(words[3])
    ####to read the input file    
    
    if  type_lag == "Yes":           
        for line in f:
            if line not in [' \n','\n','\t\n','\r\n']:
                line = line.replace("\"", "").replace("\'", "")
                words = line.split()    
                words[-1] = words[-1].strip()             
                data_chrom.append(words[0])
                start_pos.append(int(words[1]))
                end_pos.append(int(words[2]))
                types.append(words[3])
                
        f.close()
        
        
    
        for type_name in np.unique(np.array(types)):
            result_all = {}
            for jj in ['AA', 'AC', 'AT', 'AG', 'CC', 'CA', 'CT', 'CG', 'GG', 'GA', 'GC', 'GT','TT', 'TA', 'TC', 'TG']:
                result_all[jj] = 0 
                
            aa = [index for index, value in enumerate(types) if (value == type_name) ]  
            sub_data_chrom = [data_chrom[i] for i in aa] 
            sub_start_pos = [start_pos[i] for i in aa] 
            sub_end_pos = [end_pos[i] for i in aa]
            
            for chr_ID in np.unique(np.array(sub_data_chrom)):   
                aa = [index for index, value in enumerate(sub_data_chrom) if (value == chr_ID) ]  
                sub_sub_data_chrom = [sub_data_chrom[i] for i in aa] 
                sub_sub_start_pos = [sub_start_pos[i] for i in aa] 
                sub_sub_end_pos = [sub_end_pos[i] for i in aa]
                
                aa = [index for index, value in enumerate(lags) if (value == chr_ID) ] 
                if len(aa) > 0:            
                    seq = seqs[int(aa[0])]
                    for jj in range(len(sub_sub_start_pos)):
                        ref_seq = seq.seq[(sub_sub_start_pos[jj]-1):(sub_sub_end_pos[jj])]        
                        for idx in range((len(ref_seq)-1)):
                            if result_all.has_key(ref_seq[idx:(idx+2)]):
                                result_all[ref_seq[idx:(idx+2)]] = result_all[ref_seq[idx:(idx+2)]]+1                           
                            
                            
            result_all_all[type_name] = pd.DataFrame.from_dict(result_all,orient='index')       
            
    if  type_lag == "No":              
        for line in f:
            if line not in [' \n','\n','\t\n','\r\n']:
                line = line.replace("\"", "").replace("\'", "")
                words = line.split()    
                words[-1] = words[-1].strip()             
                data_chrom.append(words[0])
                start_pos.append(int(words[1]))
                end_pos.append(int(words[2]))
                
        f.close()
        

        result_all = {}
        for jj in ['AA', 'AC', 'AT', 'AG', 'CC', 'CA', 'CT', 'CG', 'GG', 'GA', 'GC', 'GT','TT', 'TA', 'TC', 'TG']:
            result_all[jj] = 0
            
        for chr_ID in np.unique(np.array(data_chrom)):   
            aa = [index for index, value in enumerate(data_chrom) if (value == chr_ID) ]  
            sub_data_chrom = [data_chrom[i] for i in aa] 
            sub_start_pos = [start_pos[i] for i in aa] 
            sub_end_pos = [end_pos[i] for i in aa]
            
            aa = [index for index, value in enumerate(lags) if (value == chr_ID) ] 
            if len(aa) > 0:            
                seq = seqs[int(aa[0])]
                for jj in range(len(sub_start_pos)):
                    ref_seq = seq.seq[(sub_start_pos[jj]-1):(sub_end_pos[jj])]        
                    for idx in range((len(ref_seq)-1)):
                        if result_all.has_key(ref_seq[idx:(idx+2)]):
                            result_all[ref_seq[idx:(idx+2)]] = result_all[ref_seq[idx:(idx+2)]]+1 
        result_all_all["Total"] = pd.DataFrame.from_dict(result_all,orient='index')  
            
            
        
if bedfile == 'NONE': 
    result_all = {}        
    for jj in ['AA', 'AC', 'AT', 'AG', 'CC', 'CA', 'CT', 'CG', 'GG', 'GA', 'GC', 'GT','TT', 'TA', 'TC', 'TG']:
        result_all[jj] = 0 
    for ii in range(len(seqs)):
        seq = seqs[ii]        
        ref_seq = seq.seq    
        for idx in range((len(ref_seq)-1)):
            if result_all.has_key(ref_seq[idx:(idx+2)]):
                result_all[ref_seq[idx:(idx+2)]] = result_all[ref_seq[idx:(idx+2)]]+1
                
                
    result_all_all["Total"] = pd.DataFrame.from_dict(result_all,orient='index')


output_name = args.outputbedFile
result_all_all = result_all_all.reset_index()

header="\t".join(result_all_all.columns).replace("index","Pair")
 
np.savetxt("%s_result.txt"%(output_name), result_all_all, fmt='%s', delimiter='\t',header=header,comments = "")  



































