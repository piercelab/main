from Bio.PDB import *
parser = PDBParser()

f = open("/Users/ragul/Desktop/tcr/data-collection/structures/align_chains/LIST.txt", "r")
for line in f:
    pdbid = str(line.strip())
    structure = parser.get_structure('TCR', line.strip())
    ppb=CaPPBuilder()

    #chain_D = structure[0]['D']    
    #for ppd in ppb.build_peptides(chain_D): 
        
       # print (">"+pdbid)
       # print ppd.get_sequence()


    chain_E = structure[0]['E']
    for ppe in ppb.build_peptides(chain_E): 
            
        print (">"+pdbid)
        print ppe.get_sequence()
    
         
f.close()
