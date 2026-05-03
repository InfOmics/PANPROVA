#!/bin/bash
echo "################################################################################"
echo "# PANPROVA                                                                     #"
echo "################################################################################"

sdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
sdir=`dirname $sdir`

if [ -z "$PANPROVA_PATH" ]; then
	echo "environment variable PANPROVA_PATH not set!"
	echo "using location: $sdir"
else
	echo "using location: $PANPROVA_PATH"
	sdir="$PANPROVA_PATH"
fi
echo "################################################################################"




psubfile="$sdir/psubmatrix.txt"
ngenomes="1000"
rseed="123456"
genevarprob="0.5"
locvarprob="0.01"
genedupprob="0.001"
gsetvarperc="0.01"
geneaddprob="0.9"
transtable="11"
# TODO: update readme to include these new parameters
genefusionprob="0.001"
subgenedupprob="0.001"
intrasubgenedupprob="0.9"
subgeneextdelprob="0.001"
mingenenumextdel="1"
maxgenenumextdel="3"
reusedeletedgenesprob="0.5"
translocationprob="0.001"
# not used
# genefissionprob="0.001"


show_usage(){
echo "Please read the README file to get a detailed description of the usage of this tool"
echo "Parameters are:"
echo "--oprefix path_to_file: output prefix"
echo "--igenomefile path_to_file: input root genome in .PEG format"
echo "[--hgtpoolfile path_to_file]: HGT pool file. If nto specified. An empty HGT pool is used. "
echo "[--psubfile path_to_file]: position substistution matrix. Default value is ${psubfile}."
echo "[--phylofile path_to_file]: file of phylogenomic tree of the generated population. If not specified a random tree is generated."
echo "[--ngenomes ]: number of genomes to be generated, if a tree is not provided. Default value is ${ngenomes}."
echo "[--rseed ]:random seed. Default value is ${rseed}."
echo "[--gene-var-prob ]: gene variaiton probability. Default value is ${genevarprob}."
echo "[--loc-var-prob ]: locus variation probability. Default value is ${locvarprob}."
echo "[--gene-dup-prob ]: gene duplciation probability. Default value is ${genedupprob}."
echo "[--gset-var-perc ]: gene set variation percentage. Default value is ${gsetvarperc}."
echo "[--gene-add-prob ]: gene add probability. Default value is ${geneaddprob}."
echo "[--tran-stable ]: translation table to be used for generating translations in GBFF files. Default value is ${transtable}."
echo "[--gene-fusion-prob ]: gene fusion probability. Default value is ${genefusionprob}."
echo "[--sub-gene-dup-prob ]: sub-gene duplication probability. Default value is ${subgenedupprob}."
echo "[--intra-sub-gene-dup-prob ]: intra-sub-gene duplication probability. Default value is ${intrasubgenedupprob}."
echo "[--sub-gene-ext-del-prob ]: extended sub-gene deletion probability. Default value is ${subgeneextdelprob}."
echo "[--min-gene-num-ext-del ]: minimum number of genes per extended deletion event. Default value is ${mingenenumextdel}."
echo "[--max-gene-num-ext-del ]: maximum number of genes per extended deletion event. Default value is ${maxgenenumextdel}."
echo "[--reuse-deleted-genes-prob ]: probability of reusing a deleted gene (vs discarding it). Default value is ${reusedeletedgenesprob}."
echo "[--translocation-prob ]: translocation probability (cut a sub-region from one gene and reinsert it into another). Default value is ${translocationprob}."
# echo "[--gene-fission-prob ]: gene fission probability. Default value is ${genefissionprob}."
}

# $@ is all command line parameters passed to the script.
# -o is for short options like -v
# -l is for long options with double dash like --version
# the comma separates different long options
# -a is for long options with single dash like -version
# ,gene-fission-prob: F:
options=$(getopt -l "help,oprefix:,igenome:,hgtpool:,psub:,phylo:,ngenomes:,rseed:,gene-var-prob:,loc-var-prob:,gene-dup-prob:,gset-var-perc:,gene-add-prob:,trans-table:,gene-fusion-prob:,sub-gene-dup-prob:,intra-sub-gene-dup-prob:,sub-gene-ext-del-prob:,min-gene-num-ext-del:,max-gene-num-ext-del:,reuse-deleted-genes-prob:,translocation-prob:" -o "ho:g:H:M:P:n:r:f:l:d:j:a:e:F:D:i:E:m:x:R:T:" -a -- "$@")
# set --:
# If no arguments follow this option, then the positional parameters are unset. Otherwise, the positional parameters 
# are set to the arguments, even if some of them begin with a ‘-’.
eval set -- "$options"
while true
do
case $1 in
-h|--help) 
    show_usage
    exit 0
    ;;
-o|--oprefix)
    shift
	oprefix=$1
	echo "oprefix " $oprefix
    ;;
-g|--igenome)
    shift
	igenomefile=$1
    ;;
-H|--hgtpool)
    shift
	hgtpoolfile=$1
    ;;
-M|--psub)
    shift
	psubfile=$1
    ;;
-P|--phylo)
    shift
	phylofile=$1
    ;;
-n|--ngenomes)
    shift
	ngenomes=$1
    ;;
-r|--rseed)
    shift
	rseed=$1
    ;;
-f|--gene-var-prob)
    shift
	genevarprob=$1
    ;;
-l|--loc-var-prob)
    shift
	locvarprob=$1
    ;;
-d|--gene-dup-prob)
    shift
	genedupprob=$1
    ;;
-j|--gset-var-perc)
    shift
	gsetvarperc=$1
    ;;
-a|--gene-add-prob)
    shift
	geneaddprob=$1
    ;;
-e|--trans-table)
    shift
	transtable=$1
    ;;
-F|--gene-fusion-prob)
    shift
    genefusionprob=$1
    ;;
-D|--sub-gene-dup-prob)
    shift
    subgenedupprob=$1
    ;;
-i|--intra-sub-gene-dup-prob)
    shift
    intrasubgenedupprob=$1
    ;;
-E|--sub-gene-ext-del-prob)
    shift
    subgeneextdelprob=$1
    ;;
-m|--min-gene-num-ext-del)
    shift
    mingenenumextdel=$1
    ;;
-x|--max-gene-num-ext-del)
    shift
    maxgenenumextdel=$1
    ;;
-R|--reuse-deleted-genes-prob)
    shift
    reusedeletedgenesprob=$1
    ;;
-T|--translocation-prob)
    shift
    translocationprob=$1
    ;;
# -F|--gene-fission-prob)
#     shift
#     genefissionprob=$1
#     ;;
--)
    shift
    break;;
esac
shift
done


if [ -z "${igenomefile}" ]
then
echo "No root genome file was provided!!! Please use option --igenome path_to_file."
show_usage
exit 1
fi

if [ -z "${oprefix}" ]
then
oprefixt=$(mktemp)
rm $oprefixt
oprefix=`basename $oprefixt`
mkdir $oprefix
oprefix="./${oprefix}/run"
echo ""
echo "No output prefix was used. ${oprefix} will be used as output prefix."
fi


if [ -z "${hgtpoolfile}" ]
then
echo ""
echo "No HGT pool file was provided. Using a blank HGT pool."
echo "See command create_hgt_pool to create an HGT pool from a set of .PEG files."
hgtpoolfile="${oprefix}.hgtpool"
touch ${hgtpoolfile}
fi

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
echo "genefusionprob ${genefusionprob}"
echo "subgenedupprob ${subgenedupprob}"
echo "intrasubgenedupprob ${intrasubgenedupprob}"
echo "subgeneextdelprob ${subgeneextdelprob}"
echo "mingenenumextdel ${mingenenumextdel}"
echo "maxgenenumextdel ${maxgenenumextdel}"
echo "reusedeletedgenesprob ${reusedeletedgenesprob}"
echo "translocationprob ${translocationprob}"
echo ""
# echo "genefissionprob ${genefissionprob}"


if [ -z "${phylofile}" ]
then
phylofile="${oprefix}.tree"
echo ""
echo "################################################################################"
echo "No phylogenomic tree was provided. A random tree for ${ngenomes} will be created in ${phylofile}."
cmd="python3 ${sdir}/generate_tree_2.py ${ngenomes} ${rseed} ${phylofile}"
echo "$cmd"
$cmd
echo "--------------------------------------------------------------------------------"
echo "Converting tre to phyloxml, file ${phylofile}.xml."
echo "An image can be find at ${phylofile}.xml.png"
cmd="python3 ${sdir}/tree2phyloxml.py  ${phylofile} ${phylofile}.xml"
echo "$cmd"
$cmd
fi


echo ""
echo "################################################################################"
echo "Everything is ready! Evolving..."
cmd="${sdir}/evolve ${igenomefile} ${hgtpoolfile} ${oprefix} ${phylofile} ${psubfile} ${genevarprob} ${locvarprob} ${genedupprob} ${gsetvarperc} ${geneaddprob} ${rseed} ${genefusionprob} ${subgenedupprob} ${intrasubgenedupprob} ${subgeneextdelprob} ${mingenenumextdel} ${maxgenenumextdel} ${reusedeletedgenesprob} ${translocationprob}" # to add ${genefissionprob}
echo "$cmd"
$cmd


echo ""
echo "################################################################################"
echo "Converting to GBFF and GFF+FASTA formats..."
mkdir -p $oprefix/ogenomes/
cmd="python3 ${sdir}/pegs2gxx.py ${oprefix}.genome_sequences ${oprefix}.gene_parents ${oprefix}.genes ${oprefix}/ogenomes/ ${transtable} "
echo "$cmd"
$cmd


echo ""
echo "################################################################################"
echo "Extracting pangenomic distributions..."
cmd="python3 ${sdir}/get_pan_distrs.py  ${oprefix}.gene_parents  ${oprefix}"
echo "$cmd"
$cmd



echo "################################################################################"
echo "# PANPROVA                                                                     #"
echo "################################################################################"
echo "Finish!"
