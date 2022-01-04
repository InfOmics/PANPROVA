import sys
from Bio import Phylo
import matplotlib.pyplot as plt

ifile = sys.argv[1]
ofile = sys.argv[2]

childs = dict()
nodes = set()

for line in open(ifile):
    c = line.strip().split(' ')
    n = int(c[0])
    p = int(c[1])
    if p not in childs:
        childs[p] = list()
    childs[p].append(n)
    nodes.add(p)
    nodes.add(n)



def write_node(off, node, indent):
    off.write(('\t'*indent)+'<clade>\n')
    off.write(('\t'*indent)+'<name>'+str(node)+'</name>\n')
    if node in childs:
        for c in childs[node]:
            write_node(off,c, indent+1)
    off.write(('\t'*indent)+'</clade>\n')


off = open(ofile, 'w')

off.write('<?xml version="1.0" encoding="UTF-8"?>\n')
off.write('<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n')
off.write('<phylogeny rooted="true">')
off.write('<name>A PANPROVA tree</name>\n')

write_node(off, -1,1)

off.write('</phylogeny>\n')
off.write(' </phyloxml>\n')

off.flush()
off.close()



tree = Phylo.read(ofile, "phyloxml")
tree.ladderize()  # Flip branches so deeper clades are displayed at top
fig = plt.figure(figsize=(20, len(nodes)/10), dpi=300)
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, do_show=False)
plt.savefig(ofile+'.png')

