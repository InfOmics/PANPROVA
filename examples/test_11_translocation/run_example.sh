#!/bin/bash

echo "================================================================================"
echo "PANPROVA example 11: gene-fusion / translocation"
echo "================================================================================"

echo "================================================================================"
echo "Evolving..."
cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --ngenomes 10 --trans-table 4 --oprefix ./example --number-of-fusion-cycles 5 --gene-fusion-prob 1.0 --sub-gene-dup-prob 0 --sub-gene-ext-del-prob 0 --translocation-prob 1.0 --inversion-prob 0"
echo $cmd
date
/usr/bin/time -f"%E %M" $cmd
date
