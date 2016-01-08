for f in /TCRmodeller/templates/fw_structures/*.pdb;  do 
b=$(basename "$f"); 
echo $b;  
~/Rosetta/main/source/bin/score.macosclangrelease -s $f -renumber_pdb -per_chain_renumbering -out:output -overwrite

done
