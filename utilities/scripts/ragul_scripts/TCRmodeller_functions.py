#!/usr/bin/python 

import os, sys, gzip
sys.path.insert(0, '/TCRmodeller/programs/PyRosetta')
import rosetta

from contextlib import contextmanager
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

#to supress print statements (from modeller and rosetta)
with suppress_stdout():
    rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/database -mute basic -mute core -mute protocols -renumber_pdb -ignore_unrecognized_res -ignore_zero_occupancy false")
import rosetta.protocols.grafting as graft


import numpy as np
from subprocess import Popen, PIPE

from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()

from Bio import AlignIO
from Bio.SeqUtils import seq1
from Bio import SearchIO
from Bio import SeqIO

def orient_templates_with_tmalign(fw_template, prot_chain, oriented_chain, matrixfile, TMalign_program):
    orient_achain = Popen([TMalign_program, prot_chain, fw_template,  "-m", matrixfile], stdout=PIPE, stderr=PIPE)
    stdout, stderr = orient_achain.communicate()

    with open(matrixfile) as inf:
        lines=inf.readlines()
        line2_words = lines[2].split()
        line3_words = lines[3].split()
        line4_words = lines[4].split()

#rotation matrix
        u11=line2_words[2]
        u12=line2_words[3]
        u13=line2_words[4]
        u21=line3_words[2]
        u22=line3_words[3]
        u23=line3_words[4]
        u31=line4_words[2]
        u32=line4_words[3]
        u33=line4_words[4]
        rot=np.array([[u11,u21,u31],
                      [u12,u22,u32],
                      [u13,u23,u33]],dtype=float)

#translation vector
        t1=line2_words[1]
        t2=line3_words[1]
        t3=line4_words[1]
        tran=np.array((t1,t2,t3), 'f')

        atemplate=parser.get_structure("FIXED", prot_chain)
        moving=Selection.unfold_entities(atemplate, "A")

        for atom in moving:
            atom.transform(rot, tran)

        io.set_structure(atemplate)
        io.save(oriented_chain)

#concatente alpha and beta template after alignment
def concatente_two_pdb_files(prot_chain_a, prot_chain_b, prot_chains_ab):

    filenames = [prot_chain_a, prot_chain_b]
    with open(prot_chains_ab, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

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

def find_template_for_cdr_loop_with_local_scoring(cdr_loop_sequence, tmpdir, cdr_template_path, cdr_seq_len, similarity_cutoff):

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
    inpcdrseq = cdr_loop_sequence
    best_score = -99999
    best_seq = ""
    seqid = ""    
    cdr_db_file = cdr_template_path + "/cdrseq" + str(cdr_seq_len) + ".fasta"
    handle = open(cdr_db_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
	score = 0
	dbcdrseq = record.seq
	seq_similarity = cdr_seq_template_similarity(inpcdrseq, dbcdrseq)
#	seq_similarity = cdr_seq_template_identity(inpcdrseq, dbcdrseq)
        if float(seq_similarity) <= float(similarity_cutoff):
            #print seq_similarity, similarity_cutoff, inpcdrseq, dbcdrseq
            for x in xrange(0, len(inpcdrseq)):
                score += PAM30[int(aa_map[inpcdrseq[x:x+1]])][int(aa_map[dbcdrseq[x:x+1]])]
            if (score > best_score):
                best_score = score
                best_seq = dbcdrseq
                seqid = record.id
    return seqid , best_seq
                



def find_template_for_cdr_loop(cdr_loop_sequence, needle_program, tmpdir, cdr_template, cdr_seq_len):
    
    cdr_sc_matrix = 'EPAM30'
    gap_open_penalty = '100'
    gap_ext_penalty = '1'

    from Bio.Emboss.Applications import NeedleCommandline
    needle_achain = NeedleCommandline(cmd=needle_program, asequence=tmpdir +'/'+ cdr_loop_sequence, bsequence=cdr_template, gapopen=gap_open_penalty, gapextend=gap_ext_penalty, outfile=tmpdir +'/'+ cdr_loop_sequence + '.sc' , datafile=cdr_sc_matrix , aformat='score' , sprotein='Y' )
    stdout, stderr = needle_achain()
    #print(stdout)
    #print(stderr)
    
    from itertools import groupby
    from operator import itemgetter
    with open(tmpdir +'/'+ cdr_loop_sequence + '.sc') as fin:
    #parse emboss score file, ignore blanklines and lines start with '#'
    #remove paranthesis from the file to easily do numeric sort   
        lines = [lines.replace(')','').replace('(','').strip().split() for lines in fin if lines.strip() and not lines.startswith("#")]

    #to avoid gaps in the alignment, choose only alignments with length equal to length of cdr sequence
    #create newlist of scores with alignment lengths equal to cdr_seq_len    
    #IMPORTANT#SORT before grouping based on sequence length (3rd column in file)
    lines.sort(key=lambda x: int(x[2]))
    scgrps = groupby(lines, itemgetter(2))

    # if no alignment equal to cdr_seq_len, then default is list of all alignments to pick best score
    newlist = lines

    for ki,val in scgrps:
        if int(ki) == int(cdr_seq_len):
            newlist = list(val)

    newlist.sort(key=lambda x: float(x[3]), reverse=True)
    best_template = newlist[0][1]
    #print "best_template =" , newlist[0][1]
    #print "best_template_score =" , newlist[0][3] 
    #print "alignment length =" , newlist[0][2] 
        
    from Bio import SeqIO
    for record in SeqIO.parse(cdr_template, "fasta"):
        if record.id == best_template:
            output_handle = open(tmpdir +'/'+ 'best_template_' + cdr_loop_sequence, "w")
            SeqIO.write(record, output_handle, "fasta")
            output_handle.close()
            
    needle_cdr = NeedleCommandline(cmd=needle_program, asequence=tmpdir +'/'+ cdr_loop_sequence, bsequence=tmpdir +'/'+ 'best_template_' + cdr_loop_sequence, gapopen=gap_open_penalty, gapextend=gap_ext_penalty, outfile=tmpdir +'/'+ cdr_loop_sequence + '.aln.txt' , datafile=cdr_sc_matrix , sprotein='Y')
    stdout, stderr = needle_cdr()
    #print(stdout)
    #print(stderr)


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

def extract_cdr_loop(cdr_loop_sequence, tmpdir, template_dir):

    cdr_ali = AlignIO.read(tmpdir +'/'+ cdr_loop_sequence + '.aln.txt', "emboss")
    cdr_code = cdr_ali[1].id[0:4].lower()
    cdr_chainid = cdr_ali[1].id[5]
    cdr_seq = cdr_ali[1].seq.ungap("-")
    cdr_pos = extract_cdr_loop_from_pdb(cdr_seq, cdr_code, cdr_chainid, tmpdir, template_dir, cdr_loop_sequence)
    return cdr_pos

def extract_cdr_loop_from_pdb(cdr_seq, pdbid, chainid, tmpdir, template_dir, tag):

    cdr_loop_sequence = tag #tag
    cdr_code = pdbid
    cdr_chainid = chainid
    cdr_seq_len = len(cdr_seq)
    gzpdbfile_path =  template_dir + '/%s/pdb%s.ent.gz' %(cdr_code[1:3], cdr_code)
    gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
    template_full_struct = parser.get_structure('TCR-CDR', gzpdbfile)
    
    #for chaininfo in template_full_struct.get_chains():
        #if chaininfo.get_id() == str(cdr_chainid):
            #io.set_structure(chaininfo)
    s = template_full_struct[0][cdr_chainid]
    #save structure
    io.set_structure(s)
    prot = tmpdir +'/'+ cdr_code + '_' + cdr_chainid + '_' + cdr_loop_sequence + '.pdb'
    io.save(prot)
    #renumber pdb file using rosetta
    protpose = rosetta.Pose()
    rosetta.core.import_pose.pose_from_pdb( protpose , prot )
    prot_renumbered = tmpdir +'/'+ cdr_code + '_' + cdr_chainid + '_' + cdr_loop_sequence + '.renum.pdb'
    protpose.dump_pdb(prot_renumbered)
    template_pdb_chain = parser.get_structure('TCR-CDR', prot_renumbered)
    i = template_pdb_chain[0]
    res_list = Selection.unfold_entities(i, 'R')
    pdb_start_pos = None
    pdb_end_pos = None
    for index1, residue in enumerate(res_list):
        found_pos = 0
        kk = 0
        if cdr_seq[0] == seq1(residue.resname):
            while kk < cdr_seq_len:
                if not index1+kk >= len(res_list):
                    if cdr_seq[kk] == seq1(res_list[index1+kk].resname):
                        found_pos = found_pos + 1
                kk = kk + 1
            if found_pos == cdr_seq_len:
                pdb_start_pos = index1+1
                pdb_end_pos = index1+cdr_seq_len
    
    if (pdb_start_pos or pdb_end_pos) is None:   
        print "<br><p><h4>ERROR: CDR template sequence does not match structure - " , cdr_code , cdr_chainid , "</h4>"
        sys.exit()

    #print "pos : ", pdb_start_pos, pdb_end_pos                                                    


    #Extract using Biopython BioPDB DICE (not working for all PDB files)
    #extract(structure, cdr_chainid, pdb_start_pos, pdb_end_pos, str(tmpdir +'/'+ cdr_loop_sequence + '.pdb'))
    #extract(structure, cdr_chainid, pdb_start_pos, pdb_end_pos, 'tmp.pdb')
    
    #Extract using Rosetta grafting protocol
    #from toolbox import cleanATOM
    #cleanATOM(cdr_code + '_' + cdr_chainid + '_' + cdr_loop_sequence + '.pdb')
    #clean_prot = tmpdir +'/'+ cdr_code + '_' + cdr_chainid + '_' + cdr_loop_sequence + '.clean.pdb'
    cdrpose = rosetta.Pose()
    cdrpose = graft.return_region( protpose, pdb_start_pos, pdb_end_pos)
    cdrpose.dump_pdb(tmpdir +'/'+ cdr_loop_sequence + '.pdb')    
    return pdb_start_pos, pdb_end_pos

def struct_align_cdr_loop(template_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program):
    
    infile = cdr_pdb + 'profit.in'
    f = open(infile,"w")             
    f.write("ATOMS N,CA,C,O\n")                             
    f.write("REFERENCE "+ template_pdb +"\n")
    f.write("MOBILE " + cdr_pdb +"\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE " + `template_begin_pos` + "-" + `template_begin_pos+1` +":" + `cdr_begin_pos` + "-" + `cdr_begin_pos+1` + "\n")
    f.write("ZONE " + `template_end_pos-1` + "-" + `template_end_pos` +":" + `cdr_end_pos-1` + "-" + `cdr_end_pos` + "\n")
    f.write("FIT\n")
    f.write("WRITE "+cdr_pdb+"fitted.pdb"+ "\n")
    f.close()                                                                                                        
    processa = Popen([profit_program, '-f', infile], stdout=PIPE, stderr=PIPE)
    stdout, stderr = processa.communicate()
    #print(stdout)
    #print(stderr)

def run_blast(blast_program, query, database, evalue_cutoff, output_format, output_file, scoring_matrix, max_num_hits):
    process = Popen([blast_program, "-query", query, "-db", database, "-evalue", evalue_cutoff, "-outfmt", output_format, "-out" , output_file, "-matrix", scoring_matrix, "-max_target_seqs", max_num_hits ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stdout)
    print(stderr)
    
def extract_chain_and_seq_from_structure(structure, chainid, out_file):

    #to remove HETATM records
    class nonHetSelect(Select):
        def accept_residue(self,residue):
            if residue.id[0] == ' ':
                return 1
            else:
                return 0

    pdbfile = parser.get_structure("PDB", structure)
    mychain = pdbfile[0][chainid]
    io.set_structure(mychain)
    io.save(out_file+'.pdb', nonHetSelect())
    '''
    #save sequence in fasta format
    f2 = open(out_file+'.fa','w')            
    f2.write(">"+out_file+"\n")
    for ppe in ppb.build_peptides(mychain):
        f2.write(str(ppe.get_sequence()))
    f2.close()
    '''


def find_CDR_pos_from_emboss_alignment(input_alignment_file, pos):

    inp_CDR1_begin_pos = pos[0]
    inp_CDR1_end_pos = pos[1]
    inp_CDR2_begin_pos = pos[2]
    inp_CDR2_end_pos = pos[3]
    inp_CDR3_begin_pos = pos[4]
    inp_CDR3_end_pos = pos[5]
    inp_HV4_begin_pos = pos[6]
    inp_HV4_end_pos = pos[7]

    align_CDR1_begin_pos = 0
    align_CDR1_end_pos = 0
    align_CDR2_begin_pos = 0
    align_CDR2_end_pos = 0
    align_CDR3_begin_pos = 0
    align_CDR3_end_pos = 0
    align_HV4_begin_pos = 0
    align_HV4_end_pos = 0

    from Bio import AlignIO
    align = AlignIO.read(input_alignment_file, "emboss")
    #print align.format("fasta")

    for record in align :
        print(" ")
        #print(aligna[0].seq)                                                                                       
        #print(aligna[1].seq)
        

    start = 0
    for idx, record in enumerate(align[0].seq):
        #Note: idx indexed from zero                                        
        if not record.find("-") == 0:
            start += 1
            if start == inp_CDR1_begin_pos:
                align_CDR1_begin_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp += 1
                out_CDR1_begin_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_CDR1_end_pos:
                align_CDR1_end_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp -= 1
                out_CDR1_end_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_CDR2_begin_pos:
                align_CDR2_begin_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp += 1
                out_CDR2_begin_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_CDR2_end_pos:
                align_CDR2_end_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp -= 1
                out_CDR2_end_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_CDR3_begin_pos:
                align_CDR3_begin_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp += 1
                out_CDR3_begin_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_CDR3_end_pos:
                align_CDR3_end_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp -= 1
                out_CDR3_end_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_HV4_begin_pos:
                align_HV4_begin_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp += 1
                out_HV4_begin_pos = len(align[1].seq[:tmp].ungap("-")) + 1
            if start == inp_HV4_end_pos:
                align_HV4_end_pos = idx
                tmp = idx
                while align[1].seq[tmp] == "-":
                    tmp -= 1
                out_HV4_end_pos = len(align[1].seq[:tmp].ungap("-")) + 1
    '''
    print "ONE", align_CDR1_begin_pos,align_CDR1_end_pos,align_CDR2_begin_pos,align_CDR2_end_pos,align_CDR3_begin_pos,align_CDR3_end_pos,align_HV4_begin_pos,align_HV4_end_pos,out_CDR1_begin_pos, out_CDR1_end_pos, out_CDR2_begin_pos, out_CDR2_end_pos, out_CDR3_begin_pos, out_CDR3_end_pos, out_HV4_begin_pos, out_HV4_end_pos
    nogap_count = 0
    for idx, record in enumerate(align[1].seq):
        print idx, align_CDR3_begin_pos
        if not record.find("-") == 0:
            nogap_count += 1
        if idx == align_CDR1_begin_pos:
            out_CDR1_begin_pos = nogap_count
        if idx == align_CDR1_end_pos:
            out_CDR1_end_pos = nogap_count
        if idx == align_CDR2_begin_pos:
            out_CDR2_begin_pos = nogap_count
        if idx == align_CDR2_end_pos:
            out_CDR2_end_pos = nogap_count
        if idx == align_CDR3_begin_pos:
            tmp = idx
            while align[1].seq[tmp] == "-":
                tmp += 1
            out_CDR3_begin_pos = len(align[1].seq[:tmp].ungap("-")) + 1
        if idx == align_CDR3_end_pos:
            out_CDR3_end_pos = nogap_count
        if idx == align_HV4_begin_pos:
            out_HV4_begin_pos = nogap_count
        if idx == align_HV4_end_pos:
            out_HV4_end_pos = nogap_count

    print "TWO", align_CDR1_begin_pos,align_CDR1_end_pos,align_CDR2_begin_pos,align_CDR2_end_pos,align_CDR3_begin_pos,align_CDR3_end_pos,align_HV4_begin_pos,align_HV4_end_pos,out_CDR1_begin_pos, out_CDR1_end_pos, out_CDR2_begin_pos, out_CDR2_end_pos, out_CDR3_begin_pos, out_CDR3_end_pos, out_HV4_begin_pos, out_HV4_end_pos
    tmp = align[1].seq.ungap("-")
    print tmp[out_CDR1_begin_pos:out_CDR1_end_pos+1].upper(), tmp[out_CDR1_begin_pos].upper(), tmp[out_CDR1_end_pos].upper()
    print tmp[out_CDR2_begin_pos:out_CDR2_end_pos+1].upper(), tmp[out_CDR2_begin_pos].upper(), tmp[out_CDR2_end_pos].upper()
    print tmp[out_CDR3_begin_pos:out_CDR3_end_pos+1].upper(), tmp[out_CDR3_begin_pos].upper(), tmp[out_CDR3_end_pos].upper()

    print align[0].seq[align_CDR2_begin_pos:align_CDR2_end_pos+1].upper(), align_CDR2_begin_pos, align_CDR2_end_pos
    print align[1].seq[out_CDR2_begin_pos:out_CDR2_end_pos+1].upper(), out_CDR2_begin_pos, out_CDR2_end_pos 
    '''

    return(align_CDR1_begin_pos,align_CDR1_end_pos,align_CDR2_begin_pos,align_CDR2_end_pos,align_CDR3_begin_pos,align_CDR3_end_pos,align_HV4_begin_pos,align_HV4_end_pos,out_CDR1_begin_pos, out_CDR1_end_pos, out_CDR2_begin_pos, out_CDR2_end_pos, out_CDR3_begin_pos, out_CDR3_end_pos, out_HV4_begin_pos, out_HV4_end_pos)


def find_variable_domain_using_regex(tcr_seq, tag):
    import re
    
    if tag == 'a':
        regex = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[YF][A-Z]{13}([A-Z]{6,11})[A-Z]{15,30}[DL][A-Z]{2,3}Y[A-Z][CW][A-Z]([A-Z]{7,16}[FW])G[A-Z]G[A-Z]{0,7}[PA]*"
    elif tag == 'b':
        regex = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[Y][A-Z]{13}([A-Z]{6,11})[A-Z]{15,40}[YLF][A-Z][CW][A-Z]([A-Z]{7,17}[F])G[A-Z]G[A-Z]{0,7}[E]*"
        
    res = re.search(regex, str(tcr_seq))
    
    if res:
    	return res.group()
    else:
        return None

