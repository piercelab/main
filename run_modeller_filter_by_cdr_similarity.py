#!/usr/bin/env python

import os
import Bio
from Bio import SeqIO
import subprocess
from subprocess import Popen, PIPE


from sys import argv
script, aseq, bseq, cutoff  = argv


cdr_db_file = "/TCRmodeller/templates/CDR/cdrseq11.fasta"
handle = open(cdr_db_file, "rU")
for record in SeqIO.parse(handle, "fasta") :
    a="'%s'" % record.id
    b= "%s" % record.seq
    os.system("python /Users/ragul/Rosetta/tools/antibody/tcr/runmodellerserver.py -alpha_seq %s -beta_seq %s -run_dir /Users/ragul/Desktop/tmp -ignore_cdr1a  -ignore_cdr2a -ignore_hv4a  -ignore_cdr1b -ignore_cdr2b -ignore_hv4b -ignore_cdr3b -similarity_cutoff %s"%(aseq,bseq,cutoff))
   
