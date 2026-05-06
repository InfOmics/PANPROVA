#!/bin/bash

echo "================================================================================"
echo "PANPROVA example 9: gene-fusion / extended deletion"
echo "================================================================================"

echo "================================================================================"
echo "Evolving..."
cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --ngenomes 10 --trans-table 4 --oprefix ./example --number-of-fusion-cycles 5 --gene-fusion-prob 1.0 --sub-gene-dup-prob 0 --sub-gene-ext-del-prob 1.0 --min-gene-num-ext-del 1 --max-gene-num-ext-del 3 --reuse-deleted-genes-prob 0.5 --translocation-prob 0 --inversion-prob 0"
echo $cmd
date
/usr/bin/time -f"%E %M" $cmd
date
