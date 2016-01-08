# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
from modeller import soap_loop

log.verbose()
env = environ()
env.edat.contact_shell = -999
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # One loop from residue 19 to 28 inclusive
        return selection(self.residue_range('90:A', '102:A'))

m = MyLoop(env,
           inimodel='2CDF.trunc.pdb',   # initial model of the target
           sequence='2cdf',                 # code of the target
           #loop_assess_methods=assess.DOPE, # assess loops with DOPE
           loop_assess_methods=(soap_loop.Scorer(),assess.DOPE)) # assess with SOAP-Loop

m.loop.starting_model= 1           # index of the first loop model
m.loop.ending_model  = 3          # index of the last loop model
m.loop.md_level = refine.very_fast  # loop refinement method

m.make()
