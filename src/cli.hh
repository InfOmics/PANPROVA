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
    double gene_variation_prob;     // argv[6]
    double locus_variation_prob;    // argv[7]
    double gene_duplication_prob;   // argv[8]
    double geneset_variation;       // argv[9]
    double geneset_variation_add;   // argv[10]
    double geneset_variation_remove;// = 1 - add
    int   rand_seed;               // argv[11]
    double gene_fusion_prob;       // argv[12]
    double sub_gene_duplication_prob; // argv[13]
    double intrasub_gene_duplication_prob; // argv[14]
    double intersub_gene_duplication_prob; // 1 - argv[14]

    double sub_gene_extended_deletion_prob; // argv[15]
    // the following parameters includes the maximum number of genes that can
    // be involved in an extended deletion event.
    // the min number must be at least 1, and the max number must be at least equal to the min number.
    // if the min and max number are both set to 1, then the extended deletion event will be equivalent
    // to the deletion of a subsequence of a two consecutive genes.
    int min_gene_number_per_extended_deletion; // argv[16]
    int max_gene_number_per_extended_deletion; // argv[17]
    double reuse_deleted_genes_prob; // argv[18]
    double discard_deleted_genes_prob; // 1- argv[18]

    // ---------------------
    // Translocation parameters
    // ---------------------
    // probability that a translocation event occurs in the current genome.
    // a translocation cuts a codon-aligned sub-region from a source gene
    // (preserving its start and stop codons) and reinserts it into a different
    // target gene (also codon-aligned, between the target's start and stop
    // codons). neither gene is removed; the source shrinks, the target grows.
    double translocation_prob; // argv[19]

    double gene_fission_prob;       // argv[xx] not used
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

    // ---------------------
    // Gene fusion parameters
    // ---------------------
    config.gene_fusion_prob = atof(argv[12]);
    
    // sub-gene duplication parameters
    config.sub_gene_duplication_prob = atof(argv[13]);
    config.intrasub_gene_duplication_prob = atof(argv[14]);
    config.intersub_gene_duplication_prob = 1 - config.intrasub_gene_duplication_prob;


    // sub-gene extended deletion parameters
    config.sub_gene_extended_deletion_prob = atof(argv[15]);
    int min_gene_number_per_extended_deletion = atoi(argv[16]);
    if (min_gene_number_per_extended_deletion < 1) {
        std::cout << "Error: min_gene_number_per_extended_deletion must be at least 1\n";
        std::cout << "Setting min_gene_number_per_extended_deletion to 1\n";
        min_gene_number_per_extended_deletion = 1;
    }
    config.min_gene_number_per_extended_deletion = min_gene_number_per_extended_deletion;
    int max_gene_number_per_extended_deletion = atoi(argv[17]);
    if (max_gene_number_per_extended_deletion < min_gene_number_per_extended_deletion) {
        std::cout << "Error: max_gene_number_per_extended_deletion must be at least equal to min_gene_number_per_extended_deletion\n";
        std::cout << "Setting max_gene_number_per_extended_deletion to min_gene_number_per_extended_deletion\n";
        max_gene_number_per_extended_deletion = min_gene_number_per_extended_deletion;
    }
    config.max_gene_number_per_extended_deletion = max_gene_number_per_extended_deletion;

    config.reuse_deleted_genes_prob = atof(argv[18]);
    config.discard_deleted_genes_prob = 1 - config.reuse_deleted_genes_prob;

    // ---------------------
    // Translocation parameters
    // ---------------------
    config.translocation_prob = atof(argv[19]);

    // ---------------------
    // Gene fission parameters
    // ---------------------
    config.gene_fission_prob = 0; // atof(argv[xx])
}


void usage(std::string cmd){
    std::cout<<"Usage: "<<cmd<<" root_genome.peg hgt_pool.hgt oprefix tree.genome_parents sub_matrix \
    GENE_VARIATION_PROB LOCUS_VARIATION_PROB GENE_DUPLICATION_PROB GENESET_VARIATION GENESET_VARIATION_ADD RAND_SEED \
    GENE_FUSION_PROB SUB_GENE_DUPLICATION_PROB INTRA_SUB_GENE_DUPLICATION_PROB \
    SUB_GENE_EXTENDED_DELETION_PROB MIN_GENE_NUMBER_PER_EXTENDED_DELETION MAX_GENE_NUMBER_PER_EXTENDED_DELETION \
    REUSE_DELETED_GENES_PROB \
    TRANSLOCATION_PROB\n";
}

#endif