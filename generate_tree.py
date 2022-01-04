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

for i in range(1,nof_genomes):
    p = random.randint(0,i-1)
    parents.append( (i,p) )

#print(parents)

off = open(ofile,'w')
for c in parents:
    off.write(str(c[0])+' '+str(c[1])+'\n')
off.flush()
off.close()

