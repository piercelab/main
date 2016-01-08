#!/usr/bin/python 

from Bio import SeqIO

def cdr_seq_template_similarity(inpcdrseq, dbcdrseq):
    sequence_length = len(inpcdrseq)
    num_similar_residues = 0
    for x in xrange(0, sequence_length):
        if inpcdrseq[x:x+1] == dbcdrseq[x:x+1]:
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'K','R','H'} and dbcdrseq[x:x+1] in  {'K','R','H'}:#Charged positive                                                                 
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'D','E'} and dbcdrseq[x:x+1] in  {'D','E'}:#Charged negative                                                                         
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'F','Y','W'} and dbcdrseq[x:x+1] in  {'F','Y','W'}:#Aromatic                                                                         
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'V','I','L'} and dbcdrseq[x:x+1] in  {'V','I','L'}:#Aliphatic                                                                        
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'N','Q'} and dbcdrseq[x:x+1] in  {'N','Q'}:#Polar (Not charged)                                                                      
            num_similar_residues += 1
        elif inpcdrseq[x:x+1] in {'S','T'} and dbcdrseq[x:x+1] in  {'S','T'}:#Small                                                                                    
            num_similar_residues += 1

    similarity = float(num_similar_residues)/sequence_length
    return similarity


def cdr_seq_template_identity(inpcdrseq, dbcdrseq):
    sequence_length = len(inpcdrseq)
    num_identical_residues = 0
    for x in xrange(0, sequence_length):
        if inpcdrseq[x:x+1] == dbcdrseq[x:x+1]:
            num_identical_residues += 1
    identity = float(num_identical_residues)/sequence_length
    return identity

from sys import argv
script,inpcdrseq  = argv
cdr_seq_len = len(inpcdrseq)	
cdr_db_file = "/TCRmodeller/templates/CDR/cdrseq" + str(cdr_seq_len) + ".fasta"
handle = open(cdr_db_file, "rU")
for record in SeqIO.parse(handle, "fasta") :
        dbcdrseq = record.seq
        seq_similarity = cdr_seq_template_similarity(inpcdrseq, dbcdrseq)
        seq_identity = cdr_seq_template_identity(inpcdrseq, dbcdrseq)
	print inpcdrseq, dbcdrseq, seq_similarity,  seq_identity


