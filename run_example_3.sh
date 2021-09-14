date
echo "converting GB to PEG"
bash gb2peg.sh

mkdir -p example/list3_output

date
echo "creating HGT pool example/list3_output/hgt_pool.log"
g++ create_hgt_pool.cpp -o create_hgt_pool &&./create_hgt_pool example/list3 example/list3_output/hgt_pool > example/list3_output/hgt_pool.log

date
echo "evolving example/list3_output/evolve.log"
g++ evolve.cpp -o evolve &&./evolve example/input/escherichia_coli_O157H7.peg example/list3_output/hgt_pool example/list3_output/out > example/list3_output/evolve.log

date
echo "extracting pangenomic distribution example/list3_output/pandistr.log"
python3 get_pandistr_from_ogenes.py example/list3_output/out.genes > example/list3_output/pandistr.log
