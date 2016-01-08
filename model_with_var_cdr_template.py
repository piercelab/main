#!/usr/bin/env python

import os
import Bio
from Bio import SeqIO
import subprocess
from subprocess import Popen, PIPE

cdr_loop_seq = "VRMDSSYKLIF"
def calc_alignment_score(cdr_loop_seq, template_loop_seq):
    
    PAM30 = ([6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2],
             [-7,8,-6,-10,-8,-2,-9,-9,-2,-5,-8,0,-4,-9,-4,-3,-6,-2,-10,-8],
             [-4,-6,8,2,-11,-3,-2,-3,0,-5,-7,-1,-9,-9,-6,0,-2,-8,-4,-8],
             [-3,-10,2,8,-14,-2,2,-3,-4,-7,-12,-4,-11,-15,-8,-4,-5,-15,-11,-8],
             [-6,-8,-11,-14,10,-14,-14,-9,-7,-6,-15,-14,-13,-13,-8,-3,-8,-15,-4,-6],
             [-4,-2,-3,-2,-14,8,1,-7,1,-8,-5,-3,-4,-13,-3,-5,-5,-13,-12,-7],
             [-2,-9,-2,2,-14,1,8,-4,-5,-5,-9,-4,-7,-14,-5,-4,-6,-17,-8,-6],
             [-2,-9,-3,-3,-9,-7,-4,6,-9,-11,-10,-7,-8,-9,-6,-2,-6,-15,-14,-5],
             [-7,-2,0,-4,-7,1,-5,-9,9,-9,-6,-6,-10,-6,-4,-6,-7,-7,-3,-6],
             [-5,-5,-5,-7,-6,-8,-5,-11,-9,8,-1,-6,-1,-2,-8,-7,-2,-14,-6,2],
             [-6,-8,-7,-12,-15,-5,-9,-10,-6,-1,7,-8,1,-3,-7,-8,-7,-6,-7,-2],
             [-7,0,-1,-4,-14,-3,-4,-7,-6,-6,-8,7,-2,-14,-6,-4,-3,-12,-9,-9],
             [-5,-4,-9,-11,-13,-4,-7,-8,-10,-1,1,-2,11,-4,-8,-5,-4,-13,-11,-1],
             [-8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8],
             [-2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-7,-6,-8,-10,8,-2,-4,-14,-13,-6],
             [0,-3,0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6],
             [-1,-6,-2,-5,-8,-5,-6,-6,-7,-2,-7,-3,-4,-9,-4,0,7,-13,-6,-3],
             [-13,-2,-8,-15,-15,-13,-17,-15,-7,-14,-6,-12,-13,-4,-14,-5,-13,13,-5,-15],
             [-8,-10,-4,-11,-4,-12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7],
             [-2,-8,-8,-8,-6,-7,-6,-5,-6,2,-2,-9,-1,-8,-6,-6,-3,-15,-7,7])
    
    aa_map = {
        'A' : '0',	
        'R' : '1',
	'N' : '2',
	'D' : '3',
	'C' : '4',
	'Q' : '5',
	'E' : '6',
	'G' : '7',
	'H' : '8',
	'I' : '9',
	'L' : '10',
	'K' : '11',
	'M' : '12',
	'F' : '13',
	'P' : '14',
	'S' : '15',
	'T' : '16',
	'W' : '17',
	'Y' : '18',
	'V' : '19'
        }
    score = 0
    for x in xrange(0, len(cdr_loop_seq)):
        score += PAM30[int(aa_map[cdr_loop_seq[x:x+1]])][int(aa_map[template_loop_seq[x:x+1]])]
    return score



cdr_db_file = "/TCRmodeller/templates/CDR/cdrseq11.fasta"
handle = open(cdr_db_file, "rU")
for record in SeqIO.parse(handle, "fasta") :
    a="'%s'" % record.id
    b= "%s" % record.seq
    os.system("python /Users/ragul/Rosetta/tools/antibody/tcr/runmodellerserver.py -alpha_seq 3VXQ_A.fa -beta_seq 3VXQ_B.fa -achain_template 3VXQ:A -bchain_template 3VXQ:B  -run_dir /Users/ragul/Desktop/tcr/modeling/loop_template_scores -ignore_cdr1a  -ignore_cdr2a -ignore_hv4a  -ignore_cdr1b -ignore_cdr2b -ignore_hv4b -ignore_cdr3b -cdr3a_template_id %s -cdr3a_template_seq %s "%(a,b))
    pam_score = calc_alignment_score(cdr_loop_seq, b)
    profitcmd = "/Users/ragul/profit/ProFitV3.1/src/profit -f profit.in > profit.out"
    os.system("%s"%(profitcmd))
    rmscmd ="grep RMS: profit.out | tail -1 | cut -d' ' -f5"
    #    os.system("%s"%(rmscmd))
    output = subprocess.check_output(rmscmd, shell=True)
    rms = str(output).rstrip('\r\n')
    with open("test.txt", "a") as myfile:
        myfile.write(str(a)+" "+str(b)+" "+str(pam_score)+" "+str(rms)+"\n")

