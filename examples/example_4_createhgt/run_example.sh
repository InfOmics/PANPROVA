#!/bin/bash

echo "================================================================================"
echo "PANPROVA example HGT scale: running time on in creasing number of input genomes/genes"
echo "================================================================================"
echo "Full HGT genome list is in hgt_list.txt"

for (( i=2; i<9; i++))
do
echo "--------------------------------------------------------------------------------"
echo $i
head -n $i hgt_list.txt > hgt_list_${i}
cmd="../../create_hgt_pool hgt_list_${i} hgt_pool.tmp"
echo "$cmd"
/usr/bin/time -f"%E %e %M" $cmd > hgt_list_${i}.log
done

echo "================================================================================"
