import sys


if len(sys.argv) != 3:
    print("Usage: python3 phyloxml2i i.phyloxml o.genome_parents")
    exit()

ifile = sys.argv[1]
ofile = sys.argv[2]


parents = [(0,-1)]
last_id = 0
def visit_phylotree(clade, parent_id):
    global last_id
    global parents
    for c in clade.clades:
        last_id += 1
        id = last_id
        parents.append( (id,parent_id) )
        visit_phylotree(c, id)
    pass


from Bio import Phylo

xml = Phylo.parse(ifile, "phyloxml")
for tree in xml:
    print(tree)
    print('-'*20)
    #visit_phylotree(tree, 0)
    #print(tree.root.clades)
    visit_phylotree(tree.root, 0)
    print('-'*40)
    break

print(parents)

off = open(ofile,'w')
for p in parents:
    off.write(str(p[0])+' '+str(p[1])+'\n')
off.flush()
off.close()