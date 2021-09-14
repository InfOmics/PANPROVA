date
echo "converting GB to PEG"
bash gb2peg.sh

mkdir -p example/list1_output


date
echo "creating HGT pool example/list1_output/hgt_pool.log"
g++ create_hgt_pool.cpp -o create_hgt_pool &&./create_hgt_pool example/list1 example/list1_output/hgt_pool > example/list1_output/hgt_pool.log

date
echo "evolving example/list1_output/evolve.log"
g++ evolve.cpp -o evolve &&./evolve example/input/mycoplasma_genitalium_G37.peg example/list1_output/hgt_pool example/list1_output/out > example/list1_output/evolve.log

date
echo "extracting pangenomic distribution example/list1_output/pandistr.log"
python3 get_pandistr_from_ogenes.py example/list1_output/out.genes > example/list1_output/pandistr.log
