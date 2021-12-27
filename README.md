# PANPROVA
## PANgenomic PROkaryotic eVolution of full Assemblies 

PANPROVA is a computational tool for simulating pangenomic evolution by evolving the complete genomic sequence of an ancestral isolate. 
In this way, the possibility of operating at the pre-assembly stage is enabled.
Gene set variations, sequence variation and horizontal acquisition from a pool of external genomes are the evolutionary features of the tool. 

----

PANPROVA is composed of 
* a core evolution process developed in C++
* a C++ procedure for retrieving a pool of non-redundant genes extracted from a give set of genomes
* scripts for converting GBK files to the input format of the evolution software
* 3 examples and related scripts for extracting evolution statistics

----

### Evolution procedure
The evolution procedure is written in the evolve.cpp file.
User parameters are defined on top of the file. They regard:
* nof of genomes to be generated
* probability of variation when ancestor gene is acquired
* probability of variating a nucleotide in a variated sequence
* probability of duplicating a gene
* percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
* probability that the variation is a gene add
* probability that the variation if a gene removal
* seed for random number generator

Once parameters are set, the file must be compiled by typing
`g++ -o evolve  evolve.cpp`

The procedure requires two input arguments
* a file in which the root genome is reported. The file must be in the PEG format describe in what follows
* a file containing the HGT pool


In what follows, a DNA sequence is a sequence over the alphabet A,C,G,T,N.


A file in the PEG format contains an initial line in which the entire DNA sequence of the isolate is reported.
The DNA sequence is followed by the genetic annotations.
Each annotation is a new line in the file which reports the starting and the ending coordinate of the annotation and the strand (1 or -1). 
The three numbers of the annotation are separated by a single space character.


The HGT pool file is a multiline DNA file. Each line reports the complete DNA sequence of a gene.

### HGT pool extraction procedure

The procedure to extract a set of non-redundant genetic sequences for be used as HGT pool is contained in the file `create_hgt_pool.cpp`.

It applies a k-mer based similarity computation for discarding similar genetic sequences. See [1] for details regarding the measure.

The procedure has three parameters
* the length pf the k-mers used for representing a sequence
* the similarity threshold to be used when comparing HGT genes to the ancestral genes
* the similarity threshold to be used when comparing HGT genes to HGT genes

The C++ procedure must be compiled by typing 
`g++ -o create_hgt_pool create_hgt_pool.cpp`

The compiled software takes as input
* a file containing the relative paths to PEG files from which genes must be extracted. The first genome of the list is used as ancestor
* the name of the file to which write the output (the retrieved non-redundant HGT pool)


### Scripts for converting GBK files to PEG files

The Python3 script gb2peg.py converts a given GBK (GenBank) file to a PEG file (described above).
GBK files must be complete GBK files, thus they must contain the genetic annotations plus the complete DNA sequence of the isolate.

The bash script gb2peg.sh converts the GBK files in the directory example/input. Such files are used for the provided examples.

### Running the examples
Three examples are provided via the scripts run_exmaple_1.sh, run_exmaple_2.sh and run_example_3.sh

The first example evolves a M. genitalum G37 genome. The list of genomes to be used for extracting the HGT pool is provided in the file example/list1.

The first example evolves a E. coli genome. The list of genomes to be used for extracting the HGT pool is provided in the file example/list2, and it equals the list of the first example.

The third example evolves a E. coli genome. The list of genomes to be used for extracting the HGT pool is provided in the file example/list3, and it includes the list of the first example plus a set of Eschirichia genomes belonging to different strains.

The output of the first test is produced in the folder list1_output.
The logging files are produced for the example:
* hgt_pool.log
* evolve.log
* pandistr.log

The first two logging files reports details regarding the computation of the HGT pool and the evolution simulation.
The pandistr.log report the gene families produced by the evolution.
Each gene family has a numerical identifier, starting from 0, and the set of genomes in which it appears ir enclosed by brackets.
the files ends with the pangenomic distribution after a line of = symbols.
The distribution reports for each number of genomes, the number of gene families present in exact which number of genomes.

Each genome is identified with a unique numerical identifier. The -1 is assigned to the ancestor. The genomic parental relationships are reported in the file `out.genome_parents`. Each line is a relation such that the left id is the genome and the right id is its parent.

Each gene family is identified by a numerical identifier. The complete set of produce genes, those contained in the ancestor plus those retrieved from the pool, is reported in the file out.genes.
Each entry is composed of two lines. The first line reports the id of the genome, the ordering of the gene and its coordinates. The second line reports the complete gene sequence.

Genetic parental relations are reported in the file `out.gene_parents`.

Complete genomic sequences are written in the file `out.genome_sequences`.
Each line reports the id of the genome followed by the DNA sequence.

The directory of the first example also contains the script `generate_charts.py` for computing statistics regarding the generated population and charts reporting them.


----

## License
PANPROVA is distributed under the MIT license. This means that it is free for both academic and commercial use. Note however that some third party components in RI require that you reference certain works in scientific publications.
You are free to link or use RI inside source code of your own program. If do so, please reference (cite) PANPROVA and this website. We appreciate bug fixes and would be happy to collaborate for improvements. 
[License](https://raw.githubusercontent.com/InfOmics/PANPROVA/master/LICENSE)

## Citation

## References

[1] Bonnici, Vincenzo, Rosalba Giugno, and Vincenzo Manca. "PanDelos: A dictionary-based method for pan-genome content discovery." BMC bioinformatics 19.15 (2018): 47-59.
