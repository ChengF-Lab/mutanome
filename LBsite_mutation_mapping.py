#!/usr/bin/env python
import re
import sys
import os
import math
from Bio import SwissProt

UniProt_LBsites={}
with open('Dataset/BioLiP/LBsites.txt') as fp:
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        UniProt_LBsites[temp[0]]=temp[1].split(",")

def is_LBsite(UniProt_ID,mutation_site):
    indicator=None
    if UniProt_ID in UniProt_LBsites:
        site = int(re.findall(r'\d+',mutation_site)[0])
        dis = [abs(int(re.findall(r'\d+',elem)[0])-site) for elem in UniProt_LBsites[UniProt_ID]]
        min_dis = min(dis)
        nearest_site = [elem for elem in UniProt_LBsites[UniProt_ID] if abs(int(re.findall(r'\d+',elem)[0])-site) == min_dis]
        indicator = [min_dis,nearest_site]
    return indicator

maf_files = [f for f in os.listdir('/Users/junfeizhao/Google Drive/TCGA-Junfei/Somatic Mutation/') if f.endswith('NonSilent.maf')]
for maf in maf_files:
    Cancer_type = maf.split('.')[1]
    fp = open('/Users/junfeizhao/Google Drive/TCGA-Junfei/Somatic Mutation/' + maf)
    fp.next()
    fp.next()
    fp.next()
    Header = fp.next()
    LBsites_mutations = set()
    #output.write(Header)
    for line in fp:
        temp = line[0:-1].split("\t")
        if temp[8]=='Missense_Mutation':
            gene = temp[0]
            Tumor_ID = temp[15][0:12]
            AA = temp[36]
            AA_pos = re.findall(r'\d+',AA)
            if len(AA_pos) > 0:
                pos = AA_pos[0]
                Protein_ID = temp[68]
                SIFT_score = temp[72]
                PolyPhen_score = temp[73]
                dis = is_LBsite(Protein_ID,AA[2:-1])
                if dis is not None and dis[0] <= 3:
                    LBsites_mutations.add('\t'.join([Tumor_ID,gene,Protein_ID,AA,','.join(dis[1]),str(dis[0])]))
    fp.close()
    output = open('LBsites_mutations/Mutation_mapping/' + Cancer_type + '_LBsite_mutations.txt',"w+")
    output.write('\t'.join(["SampleID","Gene","UniProt","AAchange","Nearest_LBsite","Distance"]) + '\n')
    output.write('\n'.join(list(LBsites_mutations)))
    output.close()
