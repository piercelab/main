#!/usr/bin/env python

from sys import argv

script, filename = argv

IN = open(filename)

outfilename = 'unb_'+filename
OUT = open(outfilename, 'w+')

for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'D':
            OUT.write(line)
        if line[21:22] == 'E':
            OUT.write(line)
OUT.write('TER\n')


IN.close()
OUT.close()
