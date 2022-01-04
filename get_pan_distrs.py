import sys

gene_parents_file = sys.argv[1]
oprefix = sys.argv[2]


genome_list = set()
gene_list = set()
family_list  =set()
families = dict()
family2genome = dict()

gene_parents = dict()
for line in open(gene_parents_file, 'r'):
    cc = line.strip().split(" ")
    ig = (int(cc[0]),int(cc[1]))
    og = (int(cc[2]),int(cc[3]))
    gene_parents[ig] = og

    gene_list.add(ig)
    gene_list.add(og)
    genome_list.add(int(cc[0]))
    genome_list.add(int(cc[2]))

genome_list.remove(-1)
gene_list.remove((-1,-1))



def get_gene_family(genomeid, geneid):
    p = gene_parents[(genomeid,geneid)]
    if p == (-1,-1):
        return "family_"+str(genomeid)+"_"+str(geneid)
    else:
        return get_gene_family(p[0],p[1])




for g in gene_list:
    f = get_gene_family(g[0],g[1])
    family_list.add(f)

    if f not in families:
        families[f] = set()
        family2genome[f] = set()
    families[f].add(g)
    family2genome[f].add(g[0])


genomes_list = sorted(genome_list)
gene_list = sorted(gene_list)
family_list = sorted(family_list)


off = open(oprefix+".gene_families",'w')
for f in family_list:
    off.write(f+" ")
    for g in sorted(families[f]):
        off.write("("+str(g[0])+","+str(g[1])+") ")
    off.write("\n")
off.flush()
off.close()


off = open(oprefix+".family_presence",'w')
off.write("# "+" genome_".join([str(i) for i in genome_list])+"\n")
for f in family_list:
    off.write(f+" ")
    for g in genome_list:
        if g in family2genome[f]:
            off.write('x ')
        else:
            off.write('- ')
    off.write("\n")
off.flush()
off.close()


pand = dict()
for f,genomes in family2genome.items():
    c = len(genomes)
    pand[c] = pand.get(c,0)+1


off = open(oprefix+".pan_distribution",'w')
off.write("#genomes #families\n")
for c in range(1, max(pand.keys())+1):
    off.write(str(c)+" "+str(pand.get(c,0))+"\n")
off.flush()
off.close()