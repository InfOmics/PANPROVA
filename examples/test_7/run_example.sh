#!/bin/bash

echo "================================================================================"
echo "PANPROVA example 7"
echo "================================================================================"

echo "================================================================================"
echo "Evolving..."
cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --hgtpool hgt_pool --ngenomes 10 --trans-table 4 --oprefix ./example"
echo $cmd
date
/usr/bin/time -f"%E %M" $cmd
date

echo "================================================================================"
echo "Extracting pangenomic distributions and gene familes..."
date
python3 ../../get_pan_distrs.py example.gene_parents ./example
date

echo "================================================================================"
echo "Preparing input for multiple sequence alignment..."
date
mkdir imsa
python3 get_imsa_files.py example.gene_families example.genes imsa/
date


echo "================================================================================"
echo "Performign multiple sequence alignments..."
date
mkdir omsa
for f in `ls imsa/*.fna`
do
echo "----------------------------------------"
b=`basename $f | sed s/\.fna/\.msf/g`
echo $f omsa/$b
 muscle -in $f -out omsa/$b -msf -quiet
done
date



