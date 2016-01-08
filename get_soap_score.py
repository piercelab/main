#!/usr/bin/env python                                                                                                              

import os, sys
from sys import argv

script, filename, start_pos, end_pos = argv
chainid = 'A'

from modeller import *
from modeller.scripts import complete_pdb
from modeller import soap_loop
#from modeller import soap_protein_od
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')
env.edat.contact_shell = -999
sp = soap_loop.Scorer()
#sp = soap_protein_od.Scorer()
mdl = complete_pdb(env, filename)
#atmsel = selection(mdl.chains[0])
atmsel = selection(mdl.residue_range(str(start_pos)+":"+chainid, str(end_pos)+":"+chainid))
score = atmsel.assess(sp)

print "GRGR", score, filename
