#!/usr/bin/env python

import os, sys
from sys import argv
script, pdbfile = argv

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/database")

p = pose_from_pdb(pdbfile)
for i in range(1, p.total_residue() + 1):
#    print i, p.pdb_info().pose2pdb(i),p.residue(i).name1()


    filename = str(p.pdb_info().number(i))+"_mutfile.txt"
    fo = open(filename, "w")
    print "Name of the file: ", fo.name
    str1 = "total 1\n"
    str2 = "1\n"
    if not p.residue(i).name1() == "A":
        mut_res = "A"         
    elif p.residue(i).name1() == "A":
        mut_res = "G"         

    str3 = p.residue(i).name1() + " " + str(i) + " " + mut_res

    fo.write( str1 )
    fo.write( str2 )
    fo.write( str3 )


    print p.residue(i).name1() + str(p.pdb_info().number(i)) + mut_res
