#!/usr/bin/python 

from Bio.PDB import *
parser = PDBParser()

f = open("LIST.txt", "r")
for line in f:
    pdbid = str(line.strip())
    structure = parser.get_structure('TCR', line.strip())
    for model in structure:
        for chain in model:
            print line.strip() , chain.id
            
f.close()
            
