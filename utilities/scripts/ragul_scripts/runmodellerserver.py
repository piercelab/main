#!/usr/bin/python
print "Content-type: text/html"
print 

import os, sys

import cgi
form = cgi.FieldStorage()

from Bio import AlignIO

import subprocess
from subprocess import Popen

import argparse
import gzip

#change these values if required
parser = argparse.ArgumentParser()
parser.add_argument('-fw_a_database', default='/TCRmodeller/db/blast/pdb_alpha', required=False )
parser.add_argument('-fw_b_database', default='/TCRmodeller/db/blast/pdb_beta', required=False )
parser.add_argument('-a_database', default='/TCRmodeller/db/blast/tcr_alpha', required=False )
parser.add_argument('-b_database', default='/TCRmodeller/db/blast/tcr_beta', required=False )
parser.add_argument('-alpha_seq', default='', required=False )
parser.add_argument('-beta_seq', default='', required=False )
parser.add_argument('-cdr_db', default='/TCRmodeller/templates/CDR/new/ALL_CDRS.fasta', required=False )
parser.add_argument('-cdr_db_path', default='/TCRmodeller/templates/CDR', required=False )
parser.add_argument('-run_dir', default='', required=False )
parser.add_argument('-native', default='', required=False )
parser.add_argument('-calc_rmsd', default='False', required=False )
parser.add_argument('-achain_template', default=None, required=False )
parser.add_argument('-bchain_template', default=None, required=False )
parser.add_argument('-ignore_cdr1a', action='store_true', default=False)
parser.add_argument('-ignore_cdr2a', action='store_true', default=False)
parser.add_argument('-ignore_cdr3a', action='store_true', default=False)
parser.add_argument('-ignore_hv4a', action='store_true', default=False)
parser.add_argument('-ignore_cdr1b', action='store_true', default=False)
parser.add_argument('-ignore_cdr2b', action='store_true', default=False)
parser.add_argument('-ignore_cdr3b', action='store_true', default=False)
parser.add_argument('-ignore_hv4b', action='store_true', default=False)
parser.add_argument('-cdr3a_template_id', default=None, required=False )
parser.add_argument('-cdr3a_template_seq', default=None, required=False )
parser.add_argument('-similarity_cutoff', default=1, required=False )
parser.add_argument('-identity_cutoff', default=1, required=False )

args = parser.parse_args()

#PLAB server:
'''
hmmscan_program = '/usr/bin/hmmscan'
blast_program = '/TCRmodeller/programs/blast/bin/blastp'
needle_program = '/usr/bin/needle'
fast_program = '/TCRmodeller/programs/fast/fast'
TMalign_program = '/TCRmodeller/programs/TMalign/TMalign'
profit_program = '/TCRmodeller/programs/ProFitV3.1/src/profit'
'''

#Ragul Macbook laptop:
hmmscan_program = '/TCRmodeller/programs/hmmer/hmmscan'
profit_program = '/Users/ragul/profit/ProFitV3.1/src/profit'
blast_program = '/usr/local/ncbi/blast/bin/blastp'
needle_program = '/usr/local/bin/needle'
fast_program = '/TCRmodeller/programs/fast/fast'
TMalign_program = '/TCRmodeller/programs/TMalign/TMalign'


#IBBR cluster:
'''
hmmscan_program = '/home/gowthamanr/programs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan'
profit_program = '/home/gowthamanr/profit/ProFitV3.1/src/profit'
blast_program = '/ibbr/ncbi-blast/bin/blastp'
needle_program = '/ibbr/EMBOSS/bin/needle'
TMalign_program = '/home/gowthamanr/programs/TM-align/TMalign'
'''

fw_a_database = args.fw_a_database
fw_b_database = args.fw_b_database
a_database = args.a_database
b_database = args.b_database
cdr_template_dir = args.cdr_db
cdr_template_path = args.cdr_db_path
similarity_cutoff = args.similarity_cutoff
identity_cutoff = args.identity_cutoff
fw_a_template = '/TCRmodeller/templates/fw_structures/'
fw_b_template = '/TCRmodeller/templates/fw_structures/'
a_template = '/TCRmodeller/PDB_RELEASE/pdb_structures'
b_template = '/TCRmodeller/PDB_RELEASE/pdb_structures'
alpha_hmm = '/TCRmodeller/db/hmm/tcr.alpha.hmm'
beta_hmm = '/TCRmodeller/db/hmm/tcr.beta.hmm'

emboss_scoring_matrix = "EBLOSUM62"
blast_scoring_matrix = 'BLOSUM62'

cwd = os.getcwd()
basepath = '/TCRmodeller/runs/'

if (args.run_dir):
    tmpdir = args.run_dir
    if not os.path.exists(tmpdir):
        print "<p>ERROR : Directory does not exists!<br>"
        sys.exit() # Exit if no sequence is entered 
else: 
#create and cd to unique directory
    import string, datetime, random
    basename = "tcr"
    random_tag = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(5))
    datetime_tag = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    uniquedirname = "_".join([basename, random_tag, datetime_tag])
    tmpdir = basepath + uniquedirname
    
    #tmpdir = '/TCRmodeller/runs'
   
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
os.chdir(tmpdir)

# getlist() returns a list containing the
# values of the fields with the given name
alphachain = form.getfirst('alphachain')
betachain = form.getfirst('betachain')

if (alphachain and betachain):
    #remove everying but characters from input string
    achainseq=''.join(i for i in alphachain if i.isalpha())
    bchainseq=''.join(i for i in betachain if i.isalpha())
elif (args.alpha_seq and args.beta_seq):
    with open (args.alpha_seq, "r") as myfilea:
        achainseq=myfilea.read().replace('\n', '')
    with open (args.beta_seq, "r") as myfileb:
        bchainseq=myfileb.read().replace('\n', '')
else:
    print "<p>ERROR: No sequences submitted<br>"
    sys.exit() # Exit if no sequence is entered  


from TCRmodeller_functions import *

#parse input equence to keep only variable domain
adomseq = find_variable_domain_using_regex(achainseq, 'a')
bdomseq = find_variable_domain_using_regex(bchainseq, 'b')

if adomseq is None:
    adomseq = achainseq
if bdomseq is None:
    bdomseq = bchainseq

querychaina = tmpdir +'/'+ 'achain.fa'
f = open(querychaina, "w")
f.write(">tcra\n")
f.write(adomseq)
f.close() 
querychainb = tmpdir +'/'+ 'bchain.fa'
f = open(querychainb, "w")
f.write(">tcrb\n")
f.write(bdomseq)
f.close() 

print "<p>The alpha variable domain sequence identified : %s<br>" % adomseq
print "<p>The beta variable domain sequence identified : %s<br>" % bdomseq


print "<p>Finding Tempates ...<br>"
#flush
sys.stdout.flush()


#CDR alpha
aa = find_CDRs(querychaina, hmmscan_program, tmpdir, 'a')
nogap_CDR1alpha_begin_pos = aa[0]
nogap_CDR1alpha_end_pos = aa[1]
nogap_CDR2alpha_begin_pos = aa[2]
nogap_CDR2alpha_end_pos = aa[3]
nogap_CDR3alpha_begin_pos = aa[4]
nogap_CDR3alpha_end_pos = aa[5]
nogap_HV4alpha_begin_pos = aa[6]
nogap_HV4alpha_end_pos = aa[7]
len_CDR1a = (aa[1] - aa[0]) + 1
len_CDR2a = (aa[3] - aa[2]) + 1
len_HV4a = (aa[7] - aa[6]) + 1
len_CDR3a = (aa[5] - aa[4]) + 1
cdr1a_seq = aa[8]
cdr2a_seq = aa[9]
hv4a_seq = aa[10]
cdr3a_seq = aa[11]

#CDR beta
bb = find_CDRs(querychainb, hmmscan_program, tmpdir , 'b')
nogap_CDR1beta_begin_pos = bb[0]
nogap_CDR1beta_end_pos = bb[1]
nogap_CDR2beta_begin_pos = bb[2]
nogap_CDR2beta_end_pos = bb[3]
nogap_CDR3beta_begin_pos = bb[4]
nogap_CDR3beta_end_pos = bb[5]
nogap_HV4beta_begin_pos = bb[6]
nogap_HV4beta_end_pos = bb[7]
len_CDR1b = (bb[1] - bb[0]) + 1
len_CDR2b = (bb[3] - bb[2]) + 1
len_HV4b = (bb[7] - bb[6]) + 1
len_CDR3b = (bb[5] - bb[4]) + 1
cdr1b_seq = bb[8]
cdr2b_seq = bb[9]
hv4b_seq = bb[10]
cdr3b_seq = bb[11]

# arguments for running BLAST
fw_outfilea = tmpdir +'/'+ 'fw_outa.xml' 
fw_outfileb = tmpdir +'/'+ 'fw_outb.xml' 
evaluecutoff='0.00001'
outputformat='5'
max_num_hits = '500'
#find framework template
run_blast(blast_program, querychaina, fw_a_database, evaluecutoff, outputformat, fw_outfilea, blast_scoring_matrix, max_num_hits)

#first check if there are hits(templates)found in blast
if os.path.isfile(fw_outfilea):
    from xml.etree.cElementTree import ElementTree, parse 
    doc = parse(fw_outfilea)
    patrn = doc.find("BlastOutput_iterations/Iteration/Iteration_message")
    if patrn is not None:
        if patrn.tag == 'Iteration_message' and patrn.text == 'No hits found':
            print "<p>No templates found for modeling"
            sys.exit() # Exit if are no blast hits
            
blast_score_list = [] #initialize list of lists to save blast scores
result_handle = open(fw_outfilea)
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        #print (str(alignment.hit_id), str(alignment.hit_def), hsp.score, hsp.bits)
        indi_list = [str(alignment.hit_def), hsp.score, hsp.bits ];
        blast_score_list.append(indi_list)
#print blast_score_list

run_blast(blast_program, querychainb, fw_b_database, evaluecutoff, outputformat, fw_outfileb, blast_scoring_matrix, max_num_hits)
best_score = 0
result_handle = open(fw_outfileb)
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        for set in blast_score_list:
                if set[0] == str(alignment.hit_def):
                    #print (str(alignment.hit_id), str(alignment.hit_def), hsp.score, hsp.bits)
                    #print (set[0], set[1], set[2])
                    total_score = set[1] + hsp.score
                    #print(total_score)
                    if total_score > best_score:
                        best_score = total_score
                        best_fw_template = set[0]
#print('best_score = ', best_score)                        
#print('best_fw_template = ', best_fw_template)                        
print "<p>The best framework template found was : ", best_fw_template

#extract pdb chains and concatenate them (for framework template)
fw_chainida = "D"
fw_atomfilea = best_fw_template
fw_pdbcodea = fw_chainida + best_fw_template
fw_aligncodea = fw_chainida + best_fw_template

fw_chainidb = "E"
fw_atomfileb = best_fw_template
fw_pdbcodeb = fw_chainidb + best_fw_template
fw_aligncodeb = fw_chainidb + best_fw_template

from Bio.PDB import PDBParser, PDBIO
parser = PDBParser(PERMISSIVE=1,QUIET=True)

from Bio.PDB import *
io = PDBIO()

extract_chain_and_seq_from_structure(fw_a_template+fw_atomfilea, fw_chainida, "fw_template_achain")
extract_chain_and_seq_from_structure(fw_b_template+fw_atomfileb, fw_chainidb, "fw_template_bchain")

#concatente
fw_template_a = tmpdir +'/'+ "fw_template_achain.pdb"
fw_template_b = tmpdir +'/'+ "fw_template_bchain.pdb"
fw_aligncode_ab = fw_chainida + '_' + fw_chainidb + '_' + best_fw_template
concatente_two_pdb_files(fw_template_a, fw_template_b, tmpdir +'/'+ fw_aligncode_ab)


# arguments for running BLAST
outfilea = tmpdir +'/'+ 'outa.xml' 
outfileb = tmpdir +'/'+ 'outb.xml' 
evaluecutoff='0.001'
outputformat='5'
max_num_hits = '1'


if args.achain_template is not None:
    pdbcodea = args.achain_template[0:4] 
    chainida = args.achain_template[5:6] 
elif args.achain_template is None:
    
    #find template for alpha chains
    run_blast(blast_program, querychaina, a_database, evaluecutoff, outputformat, outfilea, blast_scoring_matrix, max_num_hits)
    #parse blast output (chain a)
    result_handle = open(outfilea)
    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)
    E_VALUE_THRESH = 0.001
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evaluecutoff:
            #print('****Alignment Chain A****')
            #print (alignment.title, hsp.expect)
            #print alignment.hit_def[0:4]
            #print('length:', alignment.length)    
            #print(hsp.query[0:201] + '...')
            #print(hsp.match[0:201] + '...')
            #print(hsp.sbjct[0:201] + '...')
                pdbcodea = str(alignment.hit_def[0:4])
                chainida = str(alignment.hit_def[5:6])
                aligncodea = str(pdbcodea+chainida)
                atomfilea = str(pdbcodea+".pdb")
            #print("TEMPLATE CHAIN A", pdbcodea, chainida) 

print "<p>The best template found for alpha chain: ", pdbcodea+" "+chainida

if args.bchain_template is not None:
    pdbcodeb = args.bchain_template[0:4] 
    chainidb = args.bchain_template[5:6] 
elif args.bchain_template is None:

#find template for beta chains
    run_blast(blast_program, querychainb, b_database, evaluecutoff, outputformat, outfileb, blast_scoring_matrix, max_num_hits)
# parse blast output (chain b)
    result_handle = open(outfileb)
    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)
    E_VALUE_THRESH = 0.001
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evaluecutoff:
            #print('****Alignment Chain B****')
            #print (alignment.title, hsp.expect)
            #print alignment.hit_def[0:4]
                pdbcodeb = str(alignment.hit_def[0:4])
                chainidb = str(alignment.hit_def[5:6])
                aligncodeb = str(pdbcodeb+chainidb)
                atomfileb = str(pdbcodeb+".pdb")
            #print("TEMPLATE CHAIN B", pdbcodeb, chainidb) 

print "<p>The best template found for beta chain: ", pdbcodeb+" "+chainidb

#extract pdb chains and concatenate them (for alpha and beta templates)
gzpdbfile_a_path =  a_template + '/%s/pdb%s.ent.gz' %(pdbcodea[1:3], pdbcodea) 
gzpdbfile_a = gzip.open(gzpdbfile_a_path, 'rb')

extract_chain_and_seq_from_structure(gzpdbfile_a, chainida, "template_achain")
#renumber pdb file using rosetta                                              
pose = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( pose , tmpdir +'/template_achain.pdb' )
pose.dump_pdb(tmpdir +'/template_a_renumbered.pdb')
fa = open('template_achain.fa','w')
fa.write(">template_achain\n")     
fa.write(pose.sequence())                                                                                                                 
fa.close()  

gzpdbfile_b_path =  b_template + '/%s/pdb%s.ent.gz' %(pdbcodeb[1:3], pdbcodeb) 
gzpdbfile_b = gzip.open(gzpdbfile_b_path, 'rb')

extract_chain_and_seq_from_structure(gzpdbfile_b, chainidb, "template_bchain")
#renumber pdb file using rosetta                                              
pose = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( pose , tmpdir +'/template_bchain.pdb' )
pose.dump_pdb(tmpdir +'/template_b_renumbered.pdb')
fb = open('template_bchain.fa','w')
fb.write(">template_bchain\n")     
fb.write(pose.sequence())                                                                                                                 
fb.close()  

#orient individual chain templates using fw template
#aligncode_a = 'oriented_' + pdbfilea.get_id() + '_' + chainida;
#aligncode_b = 'oriented_' + pdbfileb.get_id() + '_' + chainidb;

template_achain_pdb = tmpdir +'/template_a_renumbered.pdb';
oriented_template_achain_pdb = tmpdir + '/oriented_template_a.pdb';

template_bchain_pdb = tmpdir +'/template_b_renumbered.pdb';
oriented_template_bchain_pdb = tmpdir + '/oriented_template_b.pdb';

#use TMalign for alignment
orient_templates_with_tmalign(fw_template_a, template_achain_pdb, oriented_template_achain_pdb, tmpdir +'/'+ 'amatrix.txt', TMalign_program)
orient_templates_with_tmalign(fw_template_b, template_bchain_pdb, oriented_template_bchain_pdb, tmpdir +'/'+ 'bmatrix.txt', TMalign_program)

#use 'FAST' for alignment
#orient_achain = Popen([fast_program, fw_template_a, template_achain_pdb, "-f", oriented_template_achain_pdb], stdout=PIPE, stderr=PIPE)
#stdout, stderr = orient_achain.communicate()
#orient_bchain = Popen([fast_program, fw_template_b, template_bchain_pdb, "-f", oriented_template_bchain_pdb], stdout=PIPE, stderr=PIPE)
#stdout, stderr = orient_bchain.communicate()

#create alignment file using EMBOSS
from Bio.Emboss.Applications import NeedleCommandline

needle_achain_ali_file = 'needle_achain.txt'
needle_achain = NeedleCommandline(cmd=needle_program, asequence=querychaina, bsequence=tmpdir +'/'+ 'template_achain.fa', gapopen=10, gapextend=1, outfile=tmpdir +'/'+ needle_achain_ali_file, datafile=emboss_scoring_matrix)#, aformat3=fasta)
stdout, stderr = needle_achain()

needle_bchain_ali_file = 'needle_bchain.txt'
needle_bchain = NeedleCommandline(cmd=needle_program, asequence=querychainb, bsequence=tmpdir +'/'+ 'template_bchain.fa', gapopen=10, gapextend=1, outfile=tmpdir +'/'+ needle_bchain_ali_file, datafile=emboss_scoring_matrix)#, aformat3=fasta)
stdout, stderr = needle_bchain()

aligna = AlignIO.read(needle_achain_ali_file, "emboss")
template_aa = find_CDR_pos_from_emboss_alignment(needle_achain_ali_file, aa)
#add +1 to end pos to include end residue
align_CDR1a_begin_pos = template_aa[0]
align_CDR1a_end_pos = template_aa[1]+1
align_CDR2a_begin_pos = template_aa[2]
align_CDR2a_end_pos = template_aa[3]+1
align_CDR3a_begin_pos = template_aa[4]
align_CDR3a_end_pos = template_aa[5]+1 
align_HV4a_begin_pos = template_aa[6]
align_HV4a_end_pos = template_aa[7]+1
template_CDR1a_begin_pos = template_aa[8]
template_CDR1a_end_pos = template_aa[9]
template_CDR2a_begin_pos = template_aa[10]
template_CDR2a_end_pos = template_aa[11]
template_CDR3a_begin_pos = template_aa[12]
template_CDR3a_end_pos = template_aa[13]
template_HV4a_begin_pos = template_aa[14]
template_HV4a_end_pos = template_aa[15]

alignb = AlignIO.read(needle_bchain_ali_file, "emboss")
template_bb = find_CDR_pos_from_emboss_alignment(needle_bchain_ali_file, bb)
align_CDR1b_begin_pos = template_bb[0]
align_CDR1b_end_pos = template_bb[1]+1
align_CDR2b_begin_pos = template_bb[2]
align_CDR2b_end_pos = template_bb[3]+1
align_CDR3b_begin_pos = template_bb[4]
align_CDR3b_end_pos = template_bb[5]+1
align_HV4b_begin_pos = template_bb[6]
align_HV4b_end_pos = template_bb[7]+1
template_CDR1b_begin_pos = template_bb[8]
template_CDR1b_end_pos = template_bb[9]
template_CDR2b_begin_pos = template_bb[10]
template_CDR2b_end_pos = template_bb[11]
template_CDR3b_begin_pos = template_bb[12]
template_CDR3b_end_pos = template_bb[13]
template_HV4b_begin_pos = template_bb[14]
template_HV4b_end_pos = template_bb[15]

template_pdb_pose = rosetta.core.import_pose.pose_from_pdb(oriented_template_achain_pdb)
cdr1a_template_seq = cdr1a_seq
cdr2a_template_seq = cdr2a_seq
hv4a_template_seq = hv4a_seq
cdr3a_template_seq = cdr3a_seq
if not args.ignore_cdr1a:
#Find template for CDR loops
#for cdr loops in alpha chain
    cdr_loop_sequence = 'CDR1a.fa'
    template_dir  = a_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR1a
    template_begin_pos = template_CDR1a_begin_pos
    template_end_pos = template_CDR1a_end_pos
    res_cdr1a = find_template_for_cdr_loop_with_local_scoring(cdr1a_seq, tmpdir, cdr_template_path, len_CDR1a, similarity_cutoff)
    cdr1a_code = res_cdr1a[0][0:4]
    cdr1a_chainid = res_cdr1a[0][5:6]
    cdr1a_template_seq = res_cdr1a[1]
    print "<br><p>Template found for CDR1 alpha chain : ", cdr1a_code , cdr1a_chainid
    print "<br>",cdr1a_seq
    print "<br>",res_cdr1a[1]
    cdr1a_pos = extract_cdr_loop_from_pdb(res_cdr1a[1], cdr1a_code, cdr1a_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_achain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose
    

if not args.ignore_cdr2a:
    cdr_loop_sequence = 'CDR2a.fa'
    template_dir  = a_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR2a
    template_begin_pos = template_CDR2a_begin_pos
    template_end_pos = template_CDR2a_end_pos
    res_cdr2a = find_template_for_cdr_loop_with_local_scoring(cdr2a_seq, tmpdir, cdr_template_path, len_CDR2a, similarity_cutoff)
    cdr2a_template_seq = res_cdr2a[1]
    cdr2a_code = res_cdr2a[0][0:4]
    cdr2a_chainid = res_cdr2a[0][5:6]
    print "<br><p>Template found for CDR2 alpha chain : ", cdr2a_code , cdr2a_chainid
    print "<br>",cdr2a_seq
    print "<br>",res_cdr2a[1]
    cdr2a_pos = extract_cdr_loop_from_pdb(res_cdr2a[1], cdr2a_code, cdr2a_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_achain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose

if not args.ignore_hv4a:
    cdr_loop_sequence = 'HV4a.fa'
    template_dir  = a_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_HV4a
    template_begin_pos = template_HV4a_begin_pos
    template_end_pos = template_HV4a_end_pos
    res_hv4a = find_template_for_cdr_loop_with_local_scoring(hv4a_seq, tmpdir, cdr_template_path, len_HV4a, similarity_cutoff)
    hv4a_template_seq = res_hv4a[1]
    hv4a_code = res_hv4a[0][0:4]
    hv4a_chainid = res_hv4a[0][5:6]
    print "<br><p>Template found for HV4 alpha chain : ", hv4a_code , hv4a_chainid
    print "<br>",hv4a_seq
    print "<br>",res_hv4a[1]
    hv4a_pos = extract_cdr_loop_from_pdb(res_hv4a[1], hv4a_code, hv4a_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_achain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose

if not args.ignore_cdr3a:
    cdr_loop_sequence = 'CDR3a.fa'
    template_dir  = a_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR3a
    template_begin_pos = template_CDR3a_begin_pos
    template_end_pos = template_CDR3a_end_pos
    res_cdr3a = find_template_for_cdr_loop_with_local_scoring(cdr3a_seq, tmpdir, cdr_template_path, len_CDR3a, similarity_cutoff)
    #res_cdr3a = (args.cdr3a_template_id, str(args.cdr3a_template_seq))
    cdr3a_template_seq = res_cdr3a[1]
    cdr3a_code = res_cdr3a[0][0:4]
    cdr3a_chainid = res_cdr3a[0][5:6]
    print "<br><p>Template found for CDR3 alpha chain : ", cdr3a_code , cdr3a_chainid
    print "<br>",cdr3a_seq
    print "<br>",res_cdr3a[1]
    cdr3a_pos = extract_cdr_loop_from_pdb(res_cdr3a[1], cdr3a_code, cdr3a_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_achain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose

template_pdb_pose.dump_pdb('grafted_achain.pdb')    

cdr1b_template_seq = cdr1b_seq
cdr2b_template_seq = cdr2b_seq
hv4b_template_seq = hv4b_seq
cdr3b_template_seq = cdr3b_seq
template_pdb_pose = rosetta.core.import_pose.pose_from_pdb(oriented_template_bchain_pdb)
if not args.ignore_cdr1b:
#for cdr loops in beta chain
    cdr_loop_sequence = 'CDR1b.fa'
    template_dir  = b_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR1b
    template_begin_pos = template_CDR1b_begin_pos
    template_end_pos = template_CDR1b_end_pos
    res_cdr1b = find_template_for_cdr_loop_with_local_scoring(cdr1b_seq, tmpdir, cdr_template_path, len_CDR1b, similarity_cutoff)
    cdr1b_template_seq = res_cdr1b[1]
    cdr1b_code = res_cdr1b[0][0:4]
    cdr1b_chainid = res_cdr1b[0][5:6]
    print "<br><p>Template found for CDR1 beta chain : ", cdr1b_code , cdr1b_chainid
    print "<br>",cdr1b_seq
    print "<br>",res_cdr1b[1]
    cdr1b_pos = extract_cdr_loop_from_pdb(res_cdr1b[1], cdr1b_code, cdr1b_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_bchain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose
        
if not args.ignore_cdr2b:
    cdr_loop_sequence = 'CDR2b.fa'
    template_dir  = b_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR2b
    template_begin_pos = template_CDR2b_begin_pos
    template_end_pos = template_CDR2b_end_pos
    res_cdr2b = find_template_for_cdr_loop_with_local_scoring(cdr2b_seq, tmpdir, cdr_template_path, len_CDR2b, similarity_cutoff)
    cdr2b_template_seq = res_cdr2b[1]
    cdr2b_code = res_cdr2b[0][0:4]
    cdr2b_chainid = res_cdr2b[0][5:6]
    print "<br><p>Template found for CDR2 beta chain : ", cdr2b_code , cdr2b_chainid
    print "<br>",cdr2b_seq
    print "<br>",res_cdr2b[1]
    cdr2b_pos = extract_cdr_loop_from_pdb(res_cdr2b[1], cdr2b_code, cdr2b_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_bchain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose
        
if not args.ignore_hv4b:
    cdr_loop_sequence = 'HV4b.fa'
    template_dir  = b_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_HV4b
    template_begin_pos = template_HV4b_begin_pos
    template_end_pos = template_HV4b_end_pos
    res_hv4b = find_template_for_cdr_loop_with_local_scoring(hv4b_seq, tmpdir, cdr_template_path, len_HV4b, similarity_cutoff)
    hv4b_template_seq = res_hv4b[1]
    hv4b_code = res_hv4b[0][0:4]
    hv4b_chainid = res_hv4b[0][5:6]
    print "<br><p>Template found for HV4 beta chain : ", hv4b_code , hv4b_chainid
    print "<br>",hv4b_seq
    print "<br>",res_hv4b[1]
    hv4b_pos = extract_cdr_loop_from_pdb(res_hv4b[1], hv4b_code, hv4b_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_bchain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose
        
if not args.ignore_cdr3b:
    cdr_loop_sequence = 'CDR3b.fa'
    template_dir  = b_template
    cdr_pdb = tmpdir +'/'+ cdr_loop_sequence + '.pdb'
    cdr_begin_pos = 1
    cdr_end_pos = len_CDR3b
    template_begin_pos = template_CDR3b_begin_pos
    template_end_pos = template_CDR3b_end_pos
    res_cdr3b = find_template_for_cdr_loop_with_local_scoring(cdr3b_seq, tmpdir, cdr_template_path, len_CDR3b, similarity_cutoff)
    cdr3b_template_seq = res_cdr3b[1]
    cdr3b_code = res_cdr3b[0][0:4]
    cdr3b_chainid = res_cdr3b[0][5:6]
    print "<br><p>Template found for CDR3 beta chain : ", cdr3b_code , cdr3b_chainid
    print "<br>",cdr3b_seq
    print "<br>",res_cdr3b[1]
    cdr3b_pos = extract_cdr_loop_from_pdb(res_cdr3b[1], cdr3b_code, cdr3b_chainid, tmpdir, template_dir, cdr_loop_sequence)
    struct_align_cdr_loop(oriented_template_bchain_pdb, cdr_pdb, template_begin_pos, template_end_pos, cdr_begin_pos, cdr_end_pos, tmpdir, profit_program)
    cdr_pdb_pose = rosetta.core.import_pose.pose_from_pdb(cdr_pdb+"fitted.pdb")
    graft.delete_region( template_pdb_pose, template_begin_pos, template_end_pos)
    grafted_pose = graft.insert_pose_into_pose( template_pdb_pose, cdr_pdb_pose, template_begin_pos-1, template_begin_pos)
    template_pdb_pose = grafted_pose

template_pdb_pose.dump_pdb('grafted_bchain.pdb')    

#concatente alpha and beta template after alignment
aligncode_ab = 'grafted_ab_chain'
prot_chains_ab = tmpdir +'/'+ aligncode_ab + '.pdb' 
concatente_two_pdb_files('grafted_achain.pdb', 'grafted_bchain.pdb', prot_chains_ab)    

print "<p>Running Modeller ...<br><br>"
#flush
sys.stdout.flush()



#create alignment file
f = open(tmpdir +'/'+ 'needle_tcr.ali', "w")
#sequence
'''
f.write(">P1;tcr\n")
f.write("sequence::     : :     : :::-1.00:-1.00\n")
f.write(str(aligna[0].seq)+"/"+str(alignb[0].seq)+"*\n")
#structure
f.write(">P1;"+aligncode_ab+"\n")
#f.write("structure:"+aligncode_ab+":FIRST:"+chainida+":LAST:"+chainidb+":::-1.00:-1.00\n")
f.write("structureX:"+aligncode_ab+":   FIRST:"+chainida+":LAST:"+chainidb+":::-1.00:-1.00\n")
f.write(str(aligna[1].seq)+"/"+str(alignb[1].seq)+"*\n")
'''
#test
f.write(">P1;tcr\n")
f.write("sequence::     : :     : :::-1.00:-1.00\n")
f.write(str(aligna[0].seq)[:align_CDR1a_begin_pos] + str(cdr1a_seq) + str(aligna[0].seq)[align_CDR1a_end_pos:align_CDR2a_begin_pos] +  str(cdr2a_seq) + str(aligna[0].seq)[align_CDR2a_end_pos:align_HV4a_begin_pos] + str(hv4a_seq) + str(aligna[0].seq)[align_HV4a_end_pos:align_CDR3a_begin_pos] + str(cdr3a_seq) + str(aligna[0].seq)[align_CDR3a_end_pos:] + "/" + str(alignb[0].seq)[:align_CDR1b_begin_pos] + str(cdr1b_seq) + str(alignb[0].seq)[align_CDR1b_end_pos:align_CDR2b_begin_pos] +  str(cdr2b_seq) + str(alignb[0].seq)[align_CDR2b_end_pos:align_HV4b_begin_pos] + str(hv4b_seq) + str(alignb[0].seq)[align_HV4b_end_pos:align_CDR3b_begin_pos]+ str(cdr3b_seq) + str(alignb[0].seq)[align_CDR3b_end_pos:] + "*\n")

f.write(">P1;"+aligncode_ab+"\n")
#f.write("structureX:"+aligncode_ab+":   .:.:.:.:::-1.00:-1.00\n")
f.write("structure:"+aligncode_ab+":FIRST:"+'A'+":LAST:"+'A'+":::-1.00:-1.00\n")
f.write(str(aligna[1].seq)[:align_CDR1a_begin_pos] + str(cdr1a_template_seq) + str(aligna[1].seq)[align_CDR1a_end_pos:align_CDR2a_begin_pos] +  str(cdr2a_template_seq) + str(aligna[1].seq)[align_CDR2a_end_pos:align_HV4a_begin_pos] + str(hv4a_template_seq) + str(aligna[1].seq)[align_HV4a_end_pos:align_CDR3a_begin_pos] + str(cdr3a_template_seq) + str(aligna[1].seq)[align_CDR3a_end_pos:] + "/" + str(alignb[1].seq)[:align_CDR1b_begin_pos] + str(cdr1b_template_seq) + str(alignb[1].seq)[align_CDR1b_end_pos:align_CDR2b_begin_pos] +  str(cdr2b_template_seq) + str(alignb[1].seq)[align_CDR2b_end_pos:align_HV4b_begin_pos] + str(hv4b_template_seq) + str(alignb[1].seq)[align_HV4b_end_pos:align_CDR3b_begin_pos]+ str(cdr3b_template_seq) + str(alignb[1].seq)[align_CDR3b_end_pos:]+"*\n")

#f.write(str(aligna[1].seq)[:align_CDR1alpha_begin_pos] + str(aligna[1].seq)[align_CDR1alpha_begin_pos+2:align_CDR2alpha_begin_pos] + str(aligna[1].seq)[align_CDR2alpha_begin_pos:align_CDR3alpha_begin_pos] + str(aligna[1].seq)[align_CDR3alpha_begin_pos:] + "/" + str(alignb[1].seq)[:align_CDR1beta_begin_pos] + str(alignb[1].seq)[align_CDR1beta_begin_pos:align_CDR2beta_begin_pos] + str(alignb[1].seq)[align_CDR2beta_begin_pos:align_CDR3beta_begin_pos] + str(alignb[1].seq)[align_CDR3beta_begin_pos:] + "*\n")

f.close()

from modeller import *
from modeller.automodel import *
log.none()
env = environ()
env.io.atom_files_directory = ['.', '/TCRmodeller/templates/ab_structures' ]

a = automodel(env, 
              alnfile=tmpdir+ '/'+ 'needle_tcr.ali',
              knowns=(aligncode_ab),
              sequence='tcr',
              assess_methods=(assess.DOPE, assess.GA341)
              )
a.md_level = refine.very_fast
a.starting_model = 1
a.ending_model = 1


'''
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        return selection(self.residue_range('91:', '105:'))


a = loopmodel(env, 
              alnfile=tmpdir+ '/'+ 'needle_tcr.ali',
              knowns=(aligncode_ab,cdr1a_code,cdr2a_code,cdr3a_code,cdr1b_code,cdr2b_code,cdr3b_code), 
              sequence='tcr',
              assess_methods=(assess.DOPE, assess.GA341)
              )
a.md_level = None                   # No refinement of model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 4           # Last loop model
a.loop.md_level       = refine.fast # Loop model refinement level

a.starting_model = 1
a.ending_model = 1
'''
with suppress_stdout():
    a.make()
# Get a list of all successfully built models from a.outputs
ok_models = filter(lambda x: x['failure'] is None, a.outputs)
# Rank the models by DOPE score
key = 'DOPE score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
# Get top model
m = ok_models[0]
#print "Top model: %s (DOPE score %.3f)" % (m['name'], m[key])
top_model = m['name']
print top_model
#Link to Download model
modelfilename=tmpdir +'/'+ top_model
print modelfilename
if os.path.isfile(modelfilename):
    print "<br><p>Finished.<br>"
    print '<br><a href="downloadmodel.py?post=%s">Download model</a><br><br><br><br>' % modelfilename
else:
    print "<p>Modeller Failed!<br>"
 
#flush
sys.stdout.flush()



'''
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum80
gap_open = -10
gap_extend = -0.5

from Bio import SeqIO
best_score = 0
print"hi0"
print cdr1a
for record in SeqIO.parse('/TCRmodeller/templates/CDR/single_alpha_cdrs.fasta', "fasta"):
    alignment = pairwise2.align.globalxx(cdr1a, record.seq, one_alignment_only=1, score_only = 1)
    if alignment > best_score:
        best_score = alignment
        print best_score
        best_alignment = alignment
    print "done \n"
    
print "hi1"
print best_alignment
'''

def calc_rmsd():
    reference_pdb = args.native
    mobile_pdb = modelfilename
    
    d_chain = cwd +'/'+ 'chainD.pdb'
    e_chain = cwd +'/'+ 'chainE.pdb'
    de_chain = cwd +'/'+ 'chainDE.pdb'
    
    fw_pdbfilea = parser.get_structure('fw_pdbcodea', reference_pdb)
    for fw_chaininfoa in fw_pdbfilea.get_chains():
        if fw_chaininfoa.get_id() == 'D':
            io.set_structure(fw_chaininfoa)
            fname = d_chain
            f = open(fname,'w')
            io.save(str(fname))
            f.close()
        if fw_chaininfoa.get_id() == 'E':
            io.set_structure(fw_chaininfoa)
            fname = e_chain
            f = open(fname,'w')
            io.save(str(fname))
            f.close()
#concatente                                                                                                      
    cmd = "cat "+d_chain+" "+e_chain+" > "+de_chain+""
    os.system(cmd)
    
    new_reference_pdb = "native.pdb"
    new_mobile_pdb = "model.pdb"
    
    rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/database -mute basic -mute core -mute protocols -renumber_pdb -per_chain_renumbering -ignore_unrecognized_res -ignore_zero_occupancy false")
    pose = rosetta.Pose()
    rosetta.core.import_pose.pose_from_pdb( pose , cwd +'/'+ 'chainDE.pdb' )
    pose.dump_pdb(new_reference_pdb)
    
    pose = rosetta.Pose()
    rosetta.core.import_pose.pose_from_pdb( pose , mobile_pdb  )
    pose.dump_pdb(new_mobile_pdb)
    
    infile = 'profit.in'
    f = open(infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("REFERENCE "+ new_reference_pdb +"\n")
    f.write("MOBILE " + new_mobile_pdb +"\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE A*:A*" + "\n")
    f.write("ZONE B*:B*" + "\n")
    f.write("FIT\n")
    f.write("ZONE CLEAR"+ "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "A"+`nogap_CDR1alpha_begin_pos` + "-A" + `nogap_CDR1alpha_end_pos` +":A" + `nogap_CDR1alpha_begin_pos` + "-A" + `nogap_CDR1alpha_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "A"+`nogap_CDR2alpha_begin_pos` + "-A" + `nogap_CDR2alpha_end_pos` +":A" + `nogap_CDR2alpha_begin_pos` + "-A" + `nogap_CDR2alpha_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "A"+`nogap_CDR3alpha_begin_pos` + "-A" + `nogap_CDR3alpha_end_pos` +":A" + `nogap_CDR3alpha_begin_pos` + "-A" + `nogap_CDR3alpha_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "A"+`nogap_HV4alpha_begin_pos` + "-A" + `nogap_HV4alpha_end_pos` +":A" + `nogap_HV4alpha_begin_pos` + "-A" + `nogap_HV4alpha_end_pos` + "\n")
    
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "B"+`nogap_CDR1beta_begin_pos` + "-B" + `nogap_CDR1beta_end_pos` +":B" + `nogap_CDR1beta_begin_pos` + "-B" + `nogap_CDR1beta_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "B"+`nogap_CDR2beta_begin_pos` + "-B" + `nogap_CDR2beta_end_pos` +":B" + `nogap_CDR2beta_begin_pos` + "-B" + `nogap_CDR2beta_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "B"+`nogap_CDR3beta_begin_pos` + "-B" + `nogap_CDR3beta_end_pos` +":B" + `nogap_CDR3beta_begin_pos` + "-B" + `nogap_CDR3beta_end_pos` + "\n")
    f.write("DELRZONE ALL" + "\n")
    f.write("RZONE " + "B"+`nogap_HV4beta_begin_pos` + "-B" + `nogap_HV4beta_end_pos` +":B" + `nogap_HV4beta_begin_pos` + "-B" + `nogap_HV4beta_end_pos` + "\n")

    f.close()
    processa = Popen([profit_program, '-f', infile ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = processa.communicate()
    print stdout
    print stderr
    

if (args.calc_rmsd == "True") or (args.calc_rmsd == "true"):
    if (not os.path.isfile(modelfilename)) and (args.native):
        print "<p>ERROR: No model structure avaialble for RMSD calculation<br>"
    elif not args.native:
        print "<p>ERROR: No Reference structure provided for RMSD calculation<br>"
    else:
        calc_rmsd()
        
