#!/usr/bin/env python

import re
from Bio import SwissProt

UniProt_ID_update = {}
UniProt_SQ = {}
for record in SwissProt.parse(open('/data/Data/UniProt/uniprot_sprot_human.dat')):
    UniProt_SQ[record.accessions[0]] = record.sequence
    for accession in record.accessions:
        UniProt_ID_update[accession] = record.accessions[0]

def is_correctsite(UniProt_ID,Pos,AA):
    indicator=0
    if UniProt_ID in UniProt_SQ and Pos <= len(UniProt_SQ[UniProt_ID]) and UniProt_SQ[UniProt_ID][Pos-1]==AA:
        indicator=1
    return indicator

UniProt_LBsites={}
with open('BioLiP_2019-02-6_nr.txt') as fp:
    for line in fp:
        line = line.strip()
        tmp = line.split("\t")
        UniProt_ID = tmp[17]
        Ligand = tmp[4]
        if UniProt_ID in UniProt_ID_update:
            UniProt_ID = UniProt_ID_update[UniProt_ID]
        Receptor_SQ = tmp[19]
        if UniProt_ID in UniProt_SQ:
            LBsites = tmp[8].split()
            for site in LBsites:
                pos = UniProt_SQ[UniProt_ID].find(Receptor_SQ[(int(site[1:])-1):min(int(site[1:])+4,len(Receptor_SQ))])
                if not pos == -1 and is_correctsite(UniProt_ID,pos+1,site[0]):
                    if UniProt_ID not in UniProt_LBsites:
                        UniProt_LBsites[UniProt_ID] = set()
                        UniProt_LBsites[UniProt_ID].add(site[0] + str(pos+1))
                        #UniProt_LBsites[UniProt_ID].add(site[0] + str(pos+1) + ':' + Ligand)
                    else:
                        UniProt_LBsites[UniProt_ID].add(site[0] + str(pos+1))
                        #UniProt_LBsites[UniProt_ID].add(site[0] + str(pos+1) + ':' + Ligand)

output = open("LBsites_with_Ligand.txt","w+")
for UniProt in UniProt_LBsites:
    output.write(UniProt + '\t' + ','.join(list(UniProt_LBsites[UniProt])) + '\n')
output.close()
