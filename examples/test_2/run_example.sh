#!/bin/bash

echo "================================================================================"
echo "PANPROVA example 2"
echo "================================================================================"

echo "HGT genome list is in hgt_list.txt"
echo "Extracting HGT pool..."
cmd="../../create_hgt_pool hgt_list.txt hgt_pool"
echo "$cmd"
date
/usr/bin/time -f"%E %M" $cmd
date

echo "================================================================================"
echo "Evolving..."
cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --hgtpool hgt_pool --ngenomes 1000 --trans-table 4 --oprefix ./example"
echo $cmd
date
/usr/bin/time -f"%E %M" $cmd
date
