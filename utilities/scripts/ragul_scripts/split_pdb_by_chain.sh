for f in *_0001.pdb;  do 
b=$(basename "$f"); 
echo $b;  

python split_pdb_by_chain.py $f

done
