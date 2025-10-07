import sys
import random

if len(sys.argv) != 4:
    print("Usage: python3 generate_tree.py nof_genomes random_seed ofile.genome_parents")
    exit()

nof_genomes = int(sys.argv[1])
r_seed = int(sys.argv[2])
ofile = sys.argv[3]

random.seed(r_seed)

parents = [ (0,-1) ]
childs_counts = {i:0 for i in range(nof_genomes)}
available_parents = [0]

for i in range(1,nof_genomes):
    p = random.choice(available_parents)
    childs_counts[p] += 1
    if childs_counts[p] > 1:
    	available_parents.remove(p)
    parents.append( (i,p) )
    available_parents.append(i)

#print(parents)

off = open(ofile,'w')
for c in parents:
    off.write(str(c[0])+' '+str(c[1])+'\n')
off.flush()
off.close()

