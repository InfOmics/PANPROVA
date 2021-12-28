import sys



gene2genome = dict()

for line in open(sys.argv[1]):
    cc = line.split(' ')[0].split(':')
    genome = int(cc[0])
    gene_family = int(cc[1])
    if gene_family not in gene2genome:
        gene2genome[gene_family] = set()
    gene2genome[gene_family].add(genome)

pandistr = dict()
for k,v in sorted(gene2genome.items()):
    print(k,v)
    l = len(v)
    pandistr[l] = pandistr.get(l,0) + 1

print('='*40)
for k,v in sorted(pandistr.items()):
    print(k,v)