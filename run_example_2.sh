date
echo "converting GB to PEG"
bash gb2peg.sh

mkdir -p example/list2_output

date
echo "creating HGT pool example/list2_output/hgt_pool.log"
g++ create_hgt_pool.cpp -o create_hgt_pool &&./create_hgt_pool example/list2 example/list2_output/hgt_pool > example/list2_output/hgt_pool.log

date
echo "evolving example/list2_output/evolve.log"
g++ evolve.cpp -o evolve &&./evolve example/input/escherichia_coli_O157H7.peg example/list2_output/hgt_pool example/list2_output/out > example/list2_output/evolve.log

date
echo "extracting pangenomic distribution example/list2_output/pandistr.log"
python3 get_pandistr_from_ogenes.py example/list2_output/out.genes > example/list2_output/pandistr.log
