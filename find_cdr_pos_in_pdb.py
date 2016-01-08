#!/usr/bin/env python

import os, sys
from sys import argv

import warnings
warnings.filterwarnings('ignore')

from Bio.PDB import PDBParser, CaPPBuilder
parser = PDBParser(PERMISSIVE=1,QUIET=True)
ppb = CaPPBuilder()
sys.path.append('/Users/ragul/Rosetta/tools/antibody/tcr')
from TCRmodeller_functions import *


tmpdir = os.getcwd() 
hmmscan_program = '/TCRmodeller/programs/hmmer/hmmscan'
script, filename = argv

tag = "a"
chainid = "A"
f1 = open('temp1.fa','w')            
pdbfile = parser.get_structure("PDB", filename)
mychain = pdbfile[0][chainid]
f1.write(">"+filename+"\n")
for ppe in ppb.build_peptides(mychain):
    f1.write(str(ppe.get_sequence()))
f1.write("\n")
f1.close()
aa = find_CDRs('temp1.fa', hmmscan_program, tmpdir, tag)

tag = "b"
chainid = "B"
f2 = open('temp2.fa','w')            
pdbfile = parser.get_structure("PDB", filename)
mychain = pdbfile[0][chainid]
f2.write(">"+filename+"\n")
for ppe in ppb.build_peptides(mychain):
    f2.write(str(ppe.get_sequence()))
f2.write("\n")
f2.close()
bb = find_CDRs('temp2.fa', hmmscan_program, tmpdir, tag)

print filename, aa[2], aa[3]

#pose = pose_from_pdb(filename)
#print filename, pose.pdb_info().pdb2pose('A',aa[4]), pose.pdb_info().pdb2pose('A',aa[5]), pose.pdb_info().pdb2pose('B',bb[4]), pose.pdb_info().pdb2pose('B',bb[5])
