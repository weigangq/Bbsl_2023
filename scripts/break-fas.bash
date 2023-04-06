#!/usr/bin/env bash

fas_file=$1
num=$(bioseq -n $fas_file)
for (( i=1; i<=$num; i++ )); do
    bioseq -p"order:$i" $fas_file > tmp.fas
    id_line=$(grep "^>" tmp.fas);
    rep=$(echo $id_line | cut -f1 -d' ' | tr -d ">")
    ge=$(echo $id_line | sed -E 's/^.+organism=(Borreliella) .+$/\1/')
    sp=$(echo $id_line | sed -E 's/^.+organism=Borreliella ([^ ]+)\] .+$/\1/')
    strain=$(echo $id_line | sed -E 's/^.+strain=(.+)\].+$/\1/' | tr ' ' '_' | tr '/' '_')
    id_out=${ge}_${sp}_${strain}_${rep}
#    echo $id_out
    cp tmp.fas $id_out.clean
done
exit;