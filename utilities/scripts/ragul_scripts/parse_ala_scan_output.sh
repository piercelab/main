#!/bin/bash                                                                                                     

CWD=$(pwd)

for i in {24..42} {57..72} {81..89} {107..138}; do
    echo $i | tr '\n' '\t'  >> test.xlsx
done
echo >> test.xlsx
for f in /Users/ragul/Desktop/tcr/computational_alanine_scan/Rossjohn_unique_tcrs/ala_scan/data/tcrpmhc_aho_num/*.pdb; do
    b=`basename $f .pdb`
    for i in {24..42} {57..72} {81..89} {107..138}; do
	grep "tagu $b" res.out | grep "E $i" | awk '{print $NF}' | tr '\n' '\t' >> test.xlsx
    done
    echo >> test.xlsx
done
