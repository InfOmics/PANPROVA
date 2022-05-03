import sys

from BCBio import GFF
from Bio import SeqIO
import Bio


if len(sys.argv) != 4:
    print("Usage: python3 gff2peg.py igff ifasta opeg")
    exit()

igff = sys.argv[1]
ifasta = sys.argv[2]
opeg = sys.argv[3]


in_seq_file = ifasta
in_seq_handle = open(in_seq_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()



sequence = ""
gene_loci = list()


in_file = igff
in_handle = open(in_file)
#for rec in GFF.parse(in_handle, base_dict=seq_dict):
#    print(rec)

c_start = 0
for record in GFF.parse(in_handle, base_dict=seq_dict):
    #if(record.seq[0] != '?'):
        print(record)
        print(record.seq[0:100])
        print('------')
        if c_start == 0:
            sequence = str(record.seq)
        else:
            c_start = len(sequence)
            sequence = "N" + str(record.sequence)
        for feature in record.features:
            print(feature)
            if feature.type == 'CDS' or feature.type == 'gene':
                if isinstance(feature.location , Bio.SeqFeature.FeatureLocation):
                    strand = feature.location.strand
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    print(start,end, strand)
                    if start > end:
                        print(start,end, strand)
                        t = start
                        start = end
                        end = t
                    gene_loci.append( (start,end,strand) )

in_handle.close()


sequence = sequence.upper()
gene_loci = sorted(gene_loci)


off =  open(sys.argv[3],'w')
off.write(sequence+"\n")
for locus in gene_loci:
    off.write(str(locus[0])+" "+str(locus[1])+" "+str(locus[2])+"\n")
off.flush()
off.close()