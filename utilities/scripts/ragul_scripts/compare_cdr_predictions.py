#!/usr/bin/env python

import sys, os, re
from subprocess import Popen, PIPE
import argparse
from Bio.PDB import *
from Bio import SearchIO

tmpdir = os.getcwd()
hmmscan_program = '/TCRmodeller/programs/hmmer/hmmscan'

parser = argparse.ArgumentParser()
parser.add_argument('-alpha_seq', default='', required=False )
parser.add_argument('-beta_seq', default='', required=False )
args = parser.parse_args()

aseqfile = args.alpha_seq
bseqfile = args.beta_seq

#aa = "ACDEFGHIKLMNPQRSTVWY"

#regex_alpha = "^[A-Z]+C[A-Z]([A-Z]{8,12}W)[YF][A-Z]{13}([A-Z]{6,11})[A-Z]{2}[RGVMH][FKLIYA][A-Z]{3}([A-Z]{7,9})[A-Z]{9}[DL][A-Z]{2,3}Y[A-Z][CW][A-Z]([A-Z]{7,16}[FW])G[A-Z]G"

#regex_beta = "^[A-Z]+C[A-Z]([A-Z]{8,12}W)[Y][A-Z]{13}([A-Z]{6,11})[A-Z]{2}[VILG][A-Z]{7}([A-Z]{7,9})[A-Z]{13}[YL][A-Z]+[CW][A-Z]([A-Z]{7,16}[F])G[A-Z]G"

regex_alpha = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[YF][A-Z]{13}([A-Z]{6,11})[A-Z]{15,30}[DL][A-Z]{2,3}Y[A-Z][CW][A-Z]([A-Z]{7,16}[FW])G[A-Z]G[A-Z]{0,7}[PA]*"

regex_beta = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[Y][A-Z]{13}([A-Z]{6,11})[A-Z]{15,40}[YLF][A-Z][CW][A-Z]([A-Z]{7,17}[F])G[A-Z]G[A-Z]{0,7}[E]*"



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

        #to include end residue, add +1 to end pos
        cdr1 = fragment.query.seq[CDR1_start_pos:CDR1_end_pos+1].ungap("-").upper()
        cdr2 = fragment.query.seq[CDR2_start_pos:CDR2_end_pos+1].ungap("-").upper()
        cdr3 = fragment.query.seq[CDR3_start_pos:CDR3_end_pos+1].ungap("-").upper()
        hv4 = fragment.query.seq[HV4_start_pos:HV4_end_pos+1].ungap("-").upper()

        f = open(tmpdir +'/'+ 'CDR1'+tag+'.fa', "w")
        f.write(">CDR1"+tag+"\n")
        f.write(str(cdr1))
        f.close()
        f = open(tmpdir +'/'+ 'CDR2'+tag+'.fa', "w")
        f.write(">CDR2"+tag+"\n")
        f.write(str(cdr2))
        f.close()
        f = open(tmpdir +'/'+ 'CDR3'+tag+'.fa', "w")
        f.write(">CDR3"+tag+"\n")
        f.write(str(cdr3))
        f.close()
        f = open(tmpdir +'/'+ 'HV4'+tag+'.fa', "w")
        f.write(">HV4"+tag+"\n")
        f.write(str(hv4))
        f.close()


        return nogap_CDR1_start_pos, nogap_CDR1_end_pos, nogap_CDR2_start_pos, nogap_CDR2_end_pos, nogap_CDR3_start_pos, nogap_CDR3_end_pos, nogap_HV4_start_pos, nogap_HV4_end_pos, cdr1, cdr2, hv4, cdr3


def find_CDRs_using_regex(tcr_seq, regex):
    from Bio import SeqIO
    handle = open(tcr_seq, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        #seq = record.seq[:126]#keep ony variable domain
        seq = record.seq
    res = re.search(regex, str(seq))
    if res:
        #return res.groups()[0],res.groups()[1],res.groups()[2],res.groups()[3]
	return res.group()
    else:
        return None


from Bio import SeqIO
handle = open(bseqfile, "rU")
for record in SeqIO.parse(handle, "fasta") :
    output_handle = open("temp.fasta", "w")
    SeqIO.write(record, output_handle, "fasta")
    output_handle.close()
    print "\n", record.id, record.seq
    one = find_CDRs_using_regex("temp.fasta", regex_beta)
    #two = find_CDRs("temp.fasta", hmmscan_program, tmpdir, 'a')
    if one is None:    
	print "None : " , one
    else:
	print "Res : " , one	
	
    '''
    print one[0],one[1],one[2],one[3]
    print two[8],two[9],two[10],two[11]
    if ( one[0] == two[8] and one[1] == two[9] and  one[2] == two[10] and  one[3] == two[11]):
        print "MATCH"
    else:
        print "NO MATCH"
    '''
