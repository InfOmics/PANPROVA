#!/bin/bash

echo "================================================================================"
echo "PANPROVA example: running time on increasing number of generated genomes"
echo "================================================================================"
nofs="10 20 50 100 200 500 1000 2000 5000 10000 20000 50000"

touch hgtpool
#oprefix="tmp"
#phylofile="phylofile"
hgtpoolfile="hgtpool"
psubfile="../../psubmatrix.txt"
igenomefile="../genomes/mycoplasma_genitalium_G37.peg"
ngenomes="1000"
rseed="123456"
genevarprob="0.5"
locvarprob="0.01"
genedupprob="0.001"
gsetvarperc="0.01"
geneaddprob="0.9"
transtable="4"

echo ""
echo "Current paramters are:"
echo "oprefix ${oprefix}"
echo "igenomefile ${igenomefile}"
echo "hgtpoolfile ${hgtpoolfile}"
echo "psubfile ${psubfile}"
echo "phylofile ${phylofile}"
echo "ngenomes ${ngenomes}"
echo "rseed ${rseed}"
echo "genevarprob ${genevarprob}"
echo "locvarprob ${locvarprob}"
echo "genedupprob ${genedupprob}"
echo "gsetvarperc ${gsetvarperc}"
echo "geneaddprob ${geneaddprob}"
echo "transtable ${transtable}"
echo ""


for nof in $nofs
do
echo "--------------------------------------------------------------------------------"
echo $nof
mkdir -p tmp_${nof}
ngenomes="$nof"
oprefix="tmp_${nof}/run"
phylofile="${oprefix}_phylofile"


#cmd="bash ../../PANPROVA.sh --igenome ../genomes/mycoplasma_genitalium_G37.peg --ngenomes ${nof} --trans-table 4 --oprefix tmp_${nof}/run"
#echo "$cmd"
#/usr/bin/time -f"%E %e %M" $cmd

cmd="python3 ../../generate_tree.py ${ngenomes} ${rseed} ${phylofile}"
echo "$cmd"
/usr/bin/time -f"%E %e %M" $cmd

cmd="../../evolve ${igenomefile} ${hgtpoolfile} ${oprefix} ${phylofile} ${psubfile} ${genevarprob} ${locvarprob} ${genedupprob} ${gsetvarperc} ${geneaddprob} ${rseed}"
echo "$cmd"
/usr/bin/time -f"%E %e %M" $cmd > evolve_${nof}.log

done
echo "================================================================================"
