#!/bin/bash                                                                                                                                                    

CWD=$(pwd)

while read -r xx yy zz; do
    fbn=`basename $xx .pdb`
    rm -f ${fbn}.profit.in
    echo "ATOMS N,CA,C,O" >> ${fbn}.profit.in
    echo "REFERENCE $xx" >> ${fbn}.profit.in
    echo "MOBILE ${fbn}.pdb.30.pdb" >> ${fbn}.profit.in
    echo "ZONE CLEAR" >> ${fbn}.profit.in
    echo "ZONE A*:A*" >> ${fbn}.profit.in
    echo "ZONE B*:B*" >> ${fbn}.profit.in
    echo "FIT" >> ${fbn}.profit.in
    echo "ZONE CLEAR" >> ${fbn}.profit.in
    echo "DELRZONE ALL" >> ${fbn}.profit.in
    echo "RZONE A${yy}-A${zz}:A${yy}-A${zz}" >> ${fbn}.profit.in
    /Users/ragul/profit/ProFitV3.1/src/profit -f ${fbn}.profit.in > ${fbn}.profit.out
    var1=$(grep RMS: ${fbn}.profit.out | tail -1 | cut -d' ' -f5)
    echo -e ${fbn}'\t'$var1 
done < $1
