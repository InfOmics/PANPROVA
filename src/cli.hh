#ifndef CLI_HH
#define CLI_HH

#include <string>
#include <iostream>

struct EvolveConfig {
    std::string root_genome;       // argv[1]
    std::string hgt_pool;          // argv[2]
    std::string output_prefix;           // argv[3]
    std::string itree;             // argv[4]
    std::string isubm;             // argv[5]
    float gene_variation_prob;     // argv[6]
    float locus_variation_prob;    // argv[7]
    float gene_duplication_prob;   // argv[8]
    float geneset_variation;       // argv[9]
    float geneset_variation_add;   // argv[10]
    float geneset_variation_remove;// = 1 - add
    int   rand_seed;               // argv[11]
    
};


// #define NOF_GENOMES 999
// #define GENE_VARIATION_PROB 0.5	//probability of variation when ancestor gene is aquired
// #define LOCUS_VARIATION_PROB 0.01 // 0.05  // probability of variating a nucleotide in a variated sequence
// #define GENE_DUPLICATION_PROB 0.001	//probability of duplicating a gene
// #define GENESET_VARIATION 0.01	//percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
// #define GENESET_VARIATION_ADD 0.9 //probability that the variation is a a gene add
// #define GENESET_VARIATION_REMOVE 0.1 //probability hat the variation if a gene removal
// #define RAND_SEED 123456789


void parse_args(int argc, char** argv, EvolveConfig &config){
    
    config.root_genome = argv[1];
    config.hgt_pool = argv[2];
    config.output_prefix = argv[3];
    config.itree = argv[4];
    config.isubm = argv[5];
    config.gene_variation_prob = atof(argv[6]);
    config.locus_variation_prob = atof(argv[7]);
    config.gene_duplication_prob = atof(argv[8]);
    config.geneset_variation = atof(argv[9]);
    config.geneset_variation_add = atof(argv[10]);
    config.geneset_variation_remove = 1.0 - config.geneset_variation_add;
    config.rand_seed = atoi(argv[11]);
}


void usage(std::string cmd){
    std::cout<<"Usage: "<<cmd<<" root_genome.peg hgt_pool.hgt oprefix tree.genome_parents sub_matrix GENE_VARIATION_PROB LOCUS_VARIATION_PROB GENE_DUPLICATION_PROB GENESET_VARIATION GENESET_VARIATION_ADD RAND_SEED\n";
}

#endif