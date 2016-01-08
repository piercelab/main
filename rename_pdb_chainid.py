#!/usr/bin/env python

from sys import argv

script, filename = argv

IN = open(filename)

outfilename = 'tmp_'+filename
OUT = open(outfilename, 'w+')

IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
		newline  = line[:21] + "D" + line[22:]
		OUT.write(newline)
OUT.write('TER\n')

IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'B':
		newline  = line[:21] + "E" + line[22:]
		OUT.write(newline)
OUT.write('TER\n')
'''
IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'C':
		newline  = line[:21] + "C" + line[22:]
		OUT.write(newline)
OUT.write('TER\n')

IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
		newline  = line[:21] + "D" + line[22:]
		OUT.write(newline)
OUT.write('TER\n')

IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'B':
		newline  = line[:21] + "E" + line[22:]
		OUT.write(newline)
OUT.write('TER\n')
'''

IN.close()
OUT.close()
