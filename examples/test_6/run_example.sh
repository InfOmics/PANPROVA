#!/bin/bash

echo "================================================================================"
echo "PANPROVA example 6"
echo "================================================================================"

#echo "HGT genome list is in hgt_list.txt"
#echo "Extracting HGT pool..."
#cmd="../../create_hgt_pool hgt_list.txt hgt_pool"
#echo "$cmd"
#date
#/usr/bin/time -f"%E %M" $cmd
#date

echo "================================================================================"
echo "Evolving..."
cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --hgtpool hgt_pool --trans-table 4 --oprefix ./example  -phylo example.tree"
echo $cmd
date
/usr/bin/time -f"%E %M" $cmd
date




git clone https://github.com/ghzuo/CVTree.git

cd CVTree
cmake ./
make
cd ..

mkdir gg
rm gglist
touch gglist
for l in `ls example/ogenomes/*.fna` 
do
echo $l
b=`basename $l | sed s/\.fna/\.ffn/g`
echo $l $b
cp $l gg/$b
echo $b >>gglist
done

CVTree/bin/cvtree -G gg  -g ffn  -i gglist


echo "genome_0 genome_1 genome_2 genome_3 genome_4 genome_5 genome_6 genome_7 genome_8 genome_9 genome_10 genome_11 genome_12 genome_13 genome_14 genome_15 genome_16 genome_17 genome_18 genome_19"  > dmm
sed s/genome\_[0-9]*\ //g dm/Hao.ffn.cv5.txt | sed 's/ *$//' | tail -n 20 >>dmm
Rscript hm0.R
mv Rplots.pdf similarities.pdf


Rscript trees0.R
mv Rplots.pdf tree.pdf
