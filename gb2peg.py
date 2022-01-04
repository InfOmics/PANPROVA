import sys
from Bio import SeqIO


sequence = ""
gene_loci = list()

c_start = 0
for record in SeqIO.parse(sys.argv[1], 'genbank'):
    if c_start == 0:
        sequence = str(record.seq)
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

            # loci = (start,end,strand)
            # print(loci, (loci[1]-loci[0])%3  )
            # gene_seq = sequence[loci[0]:loci[1]]
            # print( len(gene_seq)%3 )
            # print( str(Seq(gene_seq).translate(table=4)) )
            # print(feature.qualifiers['translation'])
            # #print(feature.translate(seq, cds=False))

sequence = sequence.upper()
gene_loci = sorted(gene_loci)


off =  open(sys.argv[2],'w')
off.write(sequence+"\n")
for locus in gene_loci:
    off.write(str(locus[0])+" "+str(locus[1])+" "+str(locus[2])+"\n")
off.flush()
off.close()

