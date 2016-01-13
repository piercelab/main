#!/usr/bin/env python

import os, sys
from sys import argv
import commands

script, infile = argv

pdbfile = os.path.splitext(os.path.basename(infile))[0]

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/database -no_optH true")

def run_ala_scan(infile, resfile, chainstomove):
     outfile = resfile+".out"
     commandline = '~/Rosetta/main/source/bin/rosetta_scripts.macosclangrelease -parser:protocol ~/Desktop/tmp/alascan_sample.xml -parser:view -inout:dbms:mode sqlite3 -inout:dbms:database_name rosetta_output.db3 -no_optH true -s %s -parser:script_vars pathtoresfile=%s chainstomove=%s -overwrite > %s' % (infile, resfile, chainstomove, outfile)
     res, output = commands.getstatusoutput(commandline)
     command2 = 'grep "protocols.features.DdGFeatures:  Residue " %s | awk \'{print $NF}\'' % (outfile)
     result, res_ddG = commands.getstatusoutput(command2)
     return res_ddG 

def print_resfile(start, end, p, infile, pdbfile, chainid):
     for x in xrange(int(start), int(end)+1):
          posenum = p.pdb_info().pdb2pose(chainid,x)
          if posenum == 0: 
               print "tagu", pdbfile, chainid, x, "NA", "NA", "NA"
               continue
          resfile = pdbfile+"."+str(x)+"."+chainid+".mutation.resfile"
          fo = open(resfile, "w")
          if not p.residue(posenum).name1() == "A":
               mut_res = "A"         
          elif p.residue(posenum).name1() == "A":
               mut_res = "G"         

          str1 = "NATRO\n"
          str2 = "EX 1 EX 2 EX 3\n"
          str3 = "START\n"
          str4 = str(x)+" "+chainid+" PIKAA "+mut_res+"\n"
          fo.write( str1+str2+str3+str4 )
          fo.close()
          res_ddG = run_ala_scan(infile, resfile, 2)
          print "tagu", pdbfile, chainid, x, p.residue(posenum).name1(), mut_res, res_ddG


p = pose_from_pdb(infile)
#Alpha domain
chainid = 'D'
print_resfile('24', '42', p, infile, pdbfile, chainid)#CDR1A
print_resfile('57', '72', p, infile, pdbfile, chainid)#CDR2A
print_resfile('81', '89', p, infile, pdbfile, chainid)#HV4A
print_resfile('107', '138', p, infile, pdbfile, chainid)#CDR3A
#Beta domain
chainid = 'D'
print_resfile('524', '542', p, infile, pdbfile, chainid)#CDR1B
print_resfile('557', '572', p, infile, pdbfile, chainid)#CDR2B
print_resfile('581', '589', p, infile, pdbfile, chainid)#HV4B
print_resfile('607', '638', p, infile, pdbfile, chainid)#CDR3B
