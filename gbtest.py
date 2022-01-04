import sys
from Bio import SeqIO
from Bio.Seq import Seq

sequence = ""
gene_loci = list()

rc = {'A':'T','C':'G','G':'C','T':'A','N':'N'}

c_start = 0
for record in SeqIO.parse(sys.argv[1], 'genbank'):
    if c_start == 0:
        sequence = str(record.seq)
        print(sequence[0:100])
    else:
        c_start = len(sequence)
        sequence = "N" + str(record.sequence)
    for feature in record.features:
        if feature.type == 'CDS':
            strand = feature.location.strand
            start = int(feature.location.start)
            end = int(feature.location.end)
            #print("!# ",start,end, strand, end-start, (end-start)%3)
            if start > end:
                print(start,end, strand)
                t = start
                start = end
                end = t
            gene_loci.append( (start,end,strand) )

            loci = (start,end,strand)
            print(loci, (loci[1]-loci[0])%3  )

            if strand == 1:
                gene_seq = sequence[loci[0]:loci[1]]
                print( len(gene_seq)%3 )

            #if strand == -1:
            else:
                gene_seq = sequence[loci[0]:loci[1]]
                print( len(gene_seq)%3 )
                t = list()
                for i in range(len(gene_seq)):
                    t.append(rc[gene_seq[ len(gene_seq)-1-i  ]])
                gene_seq = ''.join(t)


            print( str(Seq(gene_seq).translate(table=4)) )
            print(feature.qualifiers['translation'])
            #print(feature.translate(seq, cds=False))


sequence = sequence.upper()

for loci in gene_loci:
    print(loci, (loci[1]-loci[0])%3  )
    gene_seq = sequence[loci[0]:loci[1]]
    print( len(gene_seq)%3 )
    print( str(Seq(gene_seq).translate(table=4)) )


