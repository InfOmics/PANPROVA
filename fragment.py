import sys
from Bio import SeqIO
import math
import random

random.seed(123456)

class Fragment:
    def __init__(self):
        self.start = 0
        self.length = 0


from Bio import BiopythonParserWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonParserWarning)

    if len(sys.argv) != 6:
        print("Usage: python3 fragment.py ifile.gbff ofile.fasta nof_fragments min_fragment_length kept_genome_percentage[0...1]")
        exit()

    ifile = sys.argv[1]
    ofile = sys.argv[2]
    nof_fragments = int(sys.argv[3])
    min_fragment_length = int(sys.argv[4])
    kept_perc = float(sys.argv[5])


    genome_sequence = ""
    gene_loci = list()

    c_start = 0
    for record in SeqIO.parse(sys.argv[1], 'genbank'):
        if c_start == 0:
            genome_sequence = str(record.seq)
        else:
            c_start = len(genome_sequence)
            genome_sequence = "N" + str(record.sequence)
        for feature in record.features:
            if feature.type == 'CDS':
                strand = feature.location.strand
                start = int(feature.location.start)
                end = int(feature.location.end)
                if start > end:
                    print(start,end, strand)
                    t = start
                    start = end
                    end = t
                gene_loci.append( (start,end,strand, feature.qualifiers['locus_tag'][0]) )

    nuc_to_keep = math.floor(len(genome_sequence) * kept_perc)
    print("The input genome is "+str(len(genome_sequence))+" bp long, "+str(nuc_to_keep)+" will be kept.")

    if nof_fragments > nuc_to_keep:
        print("Too many fragments for "+str(nuc_to_keep)+" nucleotides")
        exit()

    if min_fragment_length * nof_fragments > nuc_to_keep:
        print("Number of fragments multiplied minimum fragment length is greater than "+str(nuc_to_keep))
        min_fragment_length = math.floor(nuc_to_keep / (nof_fragments*2))
        print(str(min_fragment_length)+" is the new minumum fragment length")

    
    fragments = list()
    for i in range(nof_fragments):
        fragments.append( Fragment() )
        fragments[i].length = min_fragment_length
    
    #for i in range(nuc_to_keep - (nof_fragments * min_fragment_length)):
    #    f = random.randint(0,nof_fragments-1)
    #    fragments[f].length += 1
    to_assigned = nuc_to_keep - (nof_fragments * min_fragment_length)
    while to_assigned > 0:
        f = random.randint(0,nof_fragments-1)
        n = random.randint(0, min(100, to_assigned))
        fragments[f].length += n
        to_assigned -= n

    min_distance = math.floor(( len(genome_sequence) - nuc_to_keep) / (nof_fragments+1) / 10)
    print("min distance",min_distance)

    shifts = [min_distance for i in range(nof_fragments)]
    to_assigned = len(genome_sequence) - nuc_to_keep - (min_distance*(nof_fragments+1))
    while to_assigned > 0:
        f = random.randint(0,nof_fragments-1)
        n = random.randint(0, min(100, to_assigned))
        shifts[f] += n
        to_assigned -= n


    covered_frag = [-1 for i in range(len(genome_sequence))]
    covered_pos = [-1 for i in range(len(genome_sequence))]
    fstarts = {i:0 for i in range(nof_fragments)}
    fstarts[-1] = -1
    
    off = open(ofile,'w')

    print('-'*40)
    print("fragments: id start end")
    s = 0
    for i in range(nof_fragments):
        s += shifts[i]
        fstarts[i] = s

        p1 =  s+fragments[i].start
        p2 =  s+fragments[i].start + fragments[i].length
        off.write(">fragment_"+str(i)+" "+str(p1)+" "+str(p2)+"\n")
        print("F:", i, p1, p2, sep=" ")
        off.write(genome_sequence[p1:p2]+"\n")
        for c in range(p1,p2):
            covered_frag[c] = i
            covered_pos[c] = fragments[i].start + c - s
        s +=  fragments[i].length

    off.flush()
    off.close()

    print('-'*40)
    print("genes: [Kept,Partial,Discared] locus strand fragment1 start_pos fragment2 end_pos ")
    for g in gene_loci:
        f1 = covered_frag[g[0]]
        f2 = covered_frag[g[1]]

        p1 = (f1 != -1)
        p2 = (f2 != -1)
        ctype = "D"
        if p1 and p2:
            ctype = "C"
        elif p1 or p2:
            ctype = "P"

        if f1 == -1:
            p1 = -1
        else:
            p1 = g[0] - fstarts[f1]
        if f2 == -1:
            p2 = -1
        else:
            p2 = g[1] - fstarts[f2]

        print("G:",ctype,g[3], g[2], f1, p1, f2, p2, sep=" ")










    

    