import sys
from Bio.Seq import Seq
import math

genome_sequences_file = sys.argv[1]
gene_parents_file = sys.argv[2]
genes_file = sys.argv[3]
oprefix = sys.argv[4]
trasntable = sys.argv[5]



gene_parents = dict()
for line in open(gene_parents_file, 'r'):
    cc = line.strip().split(" ")
    ig = (int(cc[0]),int(cc[1]))
    og = (int(cc[2]),int(cc[3]))
    gene_parents[ig] = og


def get_gene_family(genomeid, geneid):
    #print("@", genomeid, geneid)
    p = gene_parents[(genomeid,geneid)]
    #print(p)
    if p == (-1,-1):
        return "family_"+str(genomeid)+"_"+str(geneid)
    else:
        return get_gene_family(p[0],p[1])

#for k in gene_parents.keys():
#    print(k, get_gene_family(k[0],k[1]))


genome_lengths = dict()
for line in open(genome_sequences_file,'r'):
    cc = line.strip().split(" ")
    genome_lengths[ int(cc[0]) ] = len(cc[1])


make_new_file = set()

for line in open(genes_file,'r'):
    #cc = (line.strip().split(" "))[0].split(':')
    cc = (line.strip().split(" "))
    gene_sequence = cc[1]
    cc = cc[0].split(':')
    genome_id = int(cc[0])
    gene_id = int(cc[1])
    cc = (cc[2].replace('(','').replace(')','')).split(',')
    g_start = int(cc[0])+1
    g_end = int(cc[1])+1
    g_strand = int(cc[2])
    if g_strand==1:
        g_strand = '+'
    else:
        g_strand = '-'
    gene_family =  get_gene_family(genome_id,gene_id)


    #if (len(gene_sequence)%3) !=0:
    #    print(len(gene_sequence), line)

    if not genome_id in make_new_file:
        make_new_file.add(genome_id)

        off = open(oprefix+"genome_"+str(genome_id)+".gff", 'w')
        off.write("##gff-version 3\n")
        off.write("##sequence-region genome_"+str(genome_id)+" 1 "+str(genome_lengths[genome_id])+"\n")
        off.flush()
        off.close()

        off = open(oprefix+"genome_"+str(genome_id)+".gbff", 'w')
        off.write("LOCUS       genome_"+str(genome_id)+"                "+str(genome_lengths[genome_id])+" bp    DNA  linear   now\n")
        off.write("FEATURES             Location/Qualifiers\n")
        off.write("     source          1.."+str(genome_lengths[genome_id])+"\n")
        off.write('                     /organism="synthetic generated genome"\n')
        off.write('                     /mol_type="genomic DNA"\n')
        off.flush()
        off.close()

    off = open(oprefix+"genome_"+str(genome_id)+".gff", 'a')
    off.write("genome_"+str(genome_id)+"\tPANPROVA\tgene\t"+str(g_start)+"\t"+str(g_end)+"\t0.0\t"+g_strand+"\t0\tID="+str(gene_id)+";name="+gene_family+"\n")
    off.flush()
    off.close()

    off = open(oprefix+"genome_"+str(genome_id)+".gbff", 'a')
    if g_strand=='+':
        off.write('     gene            '+str(g_start)+'..'+str(g_end)+'\n')
    else:
        off.write('     gene            complement('+str(g_start)+'..'+str(g_end)+')\n')
    off.write('                     /gene="'+gene_family+'"\n')
    off.write('                     /locus_tag="'+str(gene_id)+'"\n')

    if g_strand=='+':
        off.write('     CDS             '+str(g_start)+'..'+str(g_end)+'\n')
    else:
        off.write('     CDS             complement('+str(g_start)+'..'+str(g_end)+')\n')
    off.write('                     /gene="'+gene_family+'"\n')
    off.write('                     /locus_tag="'+str(gene_id)+'"\n')
    off.write('                     /codon_start=1\n')
    off.write('                     /transl_table='+trasntable+'\n')
    off.write('                     /product="synthetic gene"\n')
    off.write('                     /translation="'+str(Seq(gene_sequence).translate(table=int(trasntable)))+'"\n')
    off.flush()
    off.close()

    

for line in open(genome_sequences_file,'r'):
    cc = line.strip().split(" ")

    genome_id = cc[0]
    genome_seq = cc[1]

    #print(oprefix+"genome_"+genome_id+".fna")

    off = open(oprefix+"genome_"+genome_id+".fna", 'w')
    off.write(">genome_"+genome_id+"\n")
    for i in range(0,len(genome_seq),80):
        off.write(genome_seq[i:i+80]+"\n")
    off.flush()
    off.close()


    off = open(oprefix+"genome_"+genome_id+".gbff", 'a')
    off.write("ORIGIN\n")
    for i in range(0,len(genome_seq),60):
        x = 0
        y = i+1
        while y > 0:
            x += 1
            y = math.floor(y/10)


        off.write( ' '*(9-x) )
        off.write(str(i+1)+" ")
        for j in range(i, i+60, 10):
            if j > i:
                off.write(" ")
            off.write(genome_seq[j:j+10])
        off.write("\n")
    off.write("//\n")
    off.flush()
    off.close()
