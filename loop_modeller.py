#!/usr/bin/env python

import os, sys
from sys import argv
from Bio import SearchIO
from Bio.PDB import PDBParser, CaPPBuilder
parser = PDBParser(PERMISSIVE=1,QUIET=True)
ppb = CaPPBuilder()
#sys.path.append('/Users/ragul/Rosetta/tools/antibody/tcr')
#from TCRmodeller_functions import *
from subprocess import Popen, PIPE


script, filename, tag, chainid = argv
tmpdir = os.getcwd() 
hmmscan_program = '/TCRmodeller/programs/hmmer/hmmscan'
profit_program = '/Users/ragul/profit/ProFitV3.1/src/profit'

f2 = open('temp.fa','w+')            

pdbfile = parser.get_structure("PDB", filename)
mychain = pdbfile[0][chainid]

f2.write(">"+filename+"\n")
for ppe in ppb.build_peptides(mychain):
    f2.write(str(ppe.get_sequence())+"\n")
f2.close()
    



def find_CDRs(tcr_seq, hmmscan_program, tmpdir, tag):

    CDR1_start_pos = 0 
    CDR1_end_pos = 0 
    CDR2_start_pos = 0 
    CDR2_end_pos = 0 
    CDR3_start_pos = 0 
    CDR3_end_pos = 0 
    HV4_start_pos = 0 
    HV4_end_pos = 0 

#CDR Definitions
#Reference: https://piercelab.ibbr.umd.edu/wiki/index.php/CDR_Definitions
    #CDR alpha
    if tag == 'a':
        CDR1_start = 24
        CDR1_end = 34
        CDR2_start = 49
        CDR2_end = 56
        CDR3_start = 90
        CDR3_end = 99
        HV4_start = 64
        HV4_end = 72
        tcr_hmm = '/TCRmodeller/db/hmm/tcr.alpha.hmm'
        scan_outfile = tmpdir +'/'+ 'scan_a.out'
    #CDR beta                                                                                                                                                          
    elif tag == 'b':
        CDR1_start = 24
        CDR1_end = 33
        CDR2_start = 48
        CDR2_end = 56
        CDR3_start = 93
        CDR3_end = 104
        HV4_start = 68
        HV4_end = 75
        tcr_hmm = '/TCRmodeller/db/hmm/tcr.beta.hmm'
        scan_outfile = tmpdir +'/'+ 'scan_b.out'
    
    scan_tcr = Popen([hmmscan_program, "-o", scan_outfile, tcr_hmm, tcr_seq ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = scan_tcr.communicate()
    #print(stdout)
    #print(stderr)

    for qresult in SearchIO.parse(scan_outfile, 'hmmer3-text'):
        fragment = qresult[0][0][0]
        start = fragment.hit_start
        for idx, record in enumerate(fragment.hit.seq):
            #if not fragment.hit.seq.find(".") == 0:
            if not record.find(".") == 0:
                start += 1
                if start == CDR1_start:
                    CDR1_start_pos = idx
                if start == CDR1_end:
                    CDR1_end_pos = idx
                if start == CDR2_start:
                    CDR2_start_pos = idx
                if start == CDR2_end:
                    CDR2_end_pos = idx
                if start == CDR3_start:
                    CDR3_start_pos = idx
                if start == CDR3_end:
                    CDR3_end_pos = idx
                if start == HV4_start:
                    HV4_start_pos = idx
                if start == HV4_end:
                    HV4_end_pos = idx
    
        query_pos = fragment.query_start
        for idx, record in enumerate(fragment.query.seq):
            #if not fragment.query.seq.find(".") == 0:
            if not record.find("-") == 0:
                query_pos += 1
            if idx == CDR1_start_pos:
                nogap_CDR1_start_pos = query_pos
            if idx == CDR1_end_pos:
                nogap_CDR1_end_pos = query_pos
            if idx == CDR2_start_pos:
                nogap_CDR2_start_pos = query_pos
            if idx == CDR2_end_pos:
                nogap_CDR2_end_pos = query_pos
            if idx == CDR3_start_pos:
                nogap_CDR3_start_pos = query_pos
            if idx == CDR3_end_pos:
                nogap_CDR3_end_pos = query_pos
            if idx == HV4_start_pos:
                nogap_HV4_start_pos = query_pos
            if idx == HV4_end_pos:
                nogap_HV4_end_pos = query_pos

        return nogap_CDR1_start_pos, nogap_CDR1_end_pos, nogap_CDR2_start_pos, nogap_CDR2_end_pos, nogap_CDR3_start_pos, nogap_CDR3_end_pos, nogap_HV4_start_pos, nogap_HV4_end_pos

def align_and_calc_rmsd(mobile, reference, start_pos, end_pos):
    profit_outfile = filename + '.profit.out'
    myoutput = open(profit_outfile, 'w')
    infile = 'tcr_profit.in'
    f = open(infile,"w")             
    f.write("ATOMS N,CA,C,O\n")                             
    f.write("REFERENCE "+ reference +"\n")
    f.write("MOBILE " + mobile +"\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE C*:C*" + "\n")
    f.write("ZONE D*:D*" + "\n")
    f.write("FIT\n")
    f.write("ZONE CLEAR"+ "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "C"+`start_pos` + "-C" + `end_pos` +":C" + `start_pos` + "-C" + `end_pos` + "\n")
    f.close()                                                                                                        
    processa = Popen([profit_program, '-f', infile], stdout=myoutput, stderr=PIPE)
    stdout, stderr = processa.communicate()
    #print(stdout)
    #print(stderr)


aa = find_CDRs('temp.fa', hmmscan_program, tmpdir, tag)
print filename, aa[4], aa[5] 
'''
# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
#from modeller import soap_loop

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # One loop from residue 19 to 28 inclusive
        return selection(self.residue_range(str(aa[4])+":"+chainid, str(aa[5])+":"+chainid))

m = MyLoop(env,
           inimodel=filename,   # initial model of the target
           sequence='TCR',                 # code of the target
           loop_assess_methods=assess.DOPEHR) # assess loops with DOPE
#          loop_assess_methods=soap_loop.Scorer()) # assess with SOAP-Loop

m.loop.starting_model= 1           # index of the first loop model
m.loop.ending_model  = 100           # index of the last loop model

# Very thorough VTFM optimization:
#m.library_schedule = autosched.slow
m.max_var_iterations = 300

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
m.repeat_optimization = 3
m.max_molpdf = 1e6

m.loop.md_level = refine.very_slow  # loop refinement method

m.make()

ok_models = filter(lambda x: x['failure'] is None, m.loop.outputs)
key = 'DOPE-HR score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
m = ok_models[0]
print "Top model: %s (DOPE score %.3f)" % (m['name'], m[key])
top_model = m['name']

#calc rmsd
align_and_calc_rmsd(top_model, filename, aa[4], aa[5])
'''
