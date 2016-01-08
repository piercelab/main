#!/usr/bin/env python

import sys, os, re
import gzip
from sys import argv
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/database -mute basic -mute core -mute protocols -renumber_pdb -ignore_unrecognized_res -ignore_zero_occupancy false")
import rosetta.protocols.grafting as graft

script, pdbcode, chainida, chainidb = argv

class nonHetSelect(Select):
    def accept_residue(self,residue):
        if residue.id[0] == ' ':
            return 1
        else:
            return 0

gzpdbfile_path =  "/TCRmodeller/PDB_RELEASE/pdb_structures" +  '/%s/pdb%s.ent.gz' %(pdbcode[1:3], pdbcode) 
gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
pdbfile = parser.get_structure("PDB", gzpdbfile)

mychaina = pdbfile[0][chainida]
io.set_structure(mychaina)
io.save('tmpa.pdb', nonHetSelect())
faseqa = ""
for ppe in ppb.build_peptides(mychaina):
    faseqa += str(ppe.get_sequence())

print "faseqa : ", faseqa 

mychainb = pdbfile[0][chainidb]
io.set_structure(mychainb)
io.save('tmpb.pdb', nonHetSelect())
faseqb = ""
for ppe in ppb.build_peptides(mychainb):
    faseqb += str(ppe.get_sequence())

print "faseqb : ", faseqb 

regexa = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[YF][A-Z]{13}([A-Z]{6,11})[A-Z]{15,30}[DL][A-Z]{2,3}Y[A-Z][CW][A-Z]([A-Z]{7,16}[FW])G[A-Z]G[A-Z]{0,7}[PA]*"
regexb = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[Y][A-Z]{13}([A-Z]{6,11})[A-Z]{15,40}[YLF][A-Z][CW][A-Z]([A-Z]{7,17}[F])G[A-Z]G[A-Z]{0,7}[E]*"
    
res = re.search(regexa, str(faseqa))
if res:
    print res.group(), res.start(), res.end()
else:
    print None

protposea = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( protposea , 'tmpa.pdb' )
apose = rosetta.Pose()
apose = graft.return_region( protposea, res.start()+1, res.end())
apose.dump_pdb("apose.pdb")    

res = re.search(regexb, str(faseqb))
if res:
    print res.group(), res.start(), res.end()
else:
    print None

protposeb = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( protposeb , 'tmpb.pdb' )
apose = rosetta.Pose()
apose = graft.return_region( protposeb, res.start()+1, res.end())
apose.dump_pdb("bpose.pdb")    

outfilename = pdbcode + '.trunc.pdb'
OUT = open(outfilename, 'w+')

INA = open("apose.pdb")
for line in INA.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
            newline  = line[:21] + "A" + line[22:]
            OUT.write(newline)
OUT.write('TER\n')
INA.close()

INB = open("bpose.pdb")
for line in INB.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
            newline  = line[:21] + "B" + line[22:]
            OUT.write(newline)
OUT.write('TER\n')
INB.close()
OUT.close()
