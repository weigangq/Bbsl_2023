#!/usr/bin/env bash

# standardize FASTA header for Emmanuel
# input: Excel table from "RickList BioSample-v2 Sheet2")
cat $1 | while read file species strain replicon complete topo; do
    echo -ne "working on $file ... "
    file_base=$(basename $file .fasta)
    new_file=${file_base}.clean
    new_id=">$replicon [organism=Borreliella $species] [strain=$strain] Borreliella $species strain $strain, $replicon $complete sequence"
    echo $new_id
    echo $new_id > $new_file
    grep -v ">" $file >> $new_file
    echo "done"
    echo
done

exit;
