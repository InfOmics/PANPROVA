#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <ctime>
#include <cstdlib>      // std::rand, std::srand
#include <cmath>        // std::ceil
#include <random>

#include <queue>

#include "src/lib/Locus.hh"
#include "src/lib/Genome.hh"
#include "src/utils/gene.hh"
#include "src/utils/substitution.hh"
#include "src/utils/randoms.hh"
#include "src/utils/codon.hh"
#include "src/cli.hh"
#include "src/lib/Chimera.hh"

// #define VERBOSE





int main(int argc, char** argv){

    #ifdef VERBOSE
    std::cout<<"----------------------------------------\n";
    std::cout<<"Starting evolution...\n";
    std::cout<<"verbose output enabled\n";
    std::cout<<"----------------------------------------\n";

    #endif

    std::cout<<argc<<"\n";
    if(argc!=15){
        usage(argv[0]);
        return 0;
    }

    EvolveConfig config;
    parse_args(argc, argv, config);
    
    init_random_seed(config.rand_seed);

    // std::string itree = argv[4];
    // std::string isunm = argv[5];
    // float GENE_VARIATION_PROB = atof(argv[6]);
    // float LOCUS_VARIATION_PROB = atof(argv[7]);
    // float GENE_DUPLICATION_PROB = atof(argv[8]);
    // float GENESET_VARIATION = atof(argv[9]);
    // float GENESET_VARIATION_ADD = atof(argv[10]);
    // float GENESET_VARIATION_REMOVE = 1.0 - GENESET_VARIATION_ADD;
    // int RAND_SEED = atoi(argv[11]);

/*#define NOF_GENOMES 999
#define GENE_VARIATION_PROB 0.5	//probability of variation when ancestor gene is aquired
#define LOCUS_VARIATION_PROB 0.01 // 0.05  // probability of variating a nucleotide in a variated sequence
#define GENE_DUPLICATION_PROB 0.001	//probability of duplicating a gene
#define GENESET_VARIATION 0.01	//percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
#define GENESET_VARIATION_ADD 0.9 //probability that the variation is a a gene add
#define GENESET_VARIATION_REMOVE 0.1 //probability hat the variation if a gene removal
#define RAND_SEED 123456789*/

    


    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading root genome...\n";
    Genome *root_genome = Genome::read_from_file(config.root_genome);

    std::cout<<"Root genome is "<<root_genome->sequence.size()<<" nucleotides long, with "<<root_genome->loci.size()<<" genetic loci\n";

    int min_root_gene_length=-1, max_root_gene_length=-1;
    for(Locus &l : root_genome->loci){
        int length = l.end - l.start;
        if(min_root_gene_length==-1 || min_root_gene_length>length){
            min_root_gene_length = length;
        }
        if(max_root_gene_length==-1 || max_root_gene_length<length){
            max_root_gene_length = length;
        }
        //std::cout<<"@@@ "<<l.start<<" "<<l.end<<" "<<length<<" "<<(length %3)<<"\n";
    }
    std::cout<<"gene lenghts are between "<<min_root_gene_length<<" and "<<max_root_gene_length<<" \n";

    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading HGT pool...\n";
    std::vector<std::string> hgt_pool;
    std::ifstream file(config.hgt_pool);
    std::string line;
    std::string gseq;
    while(std::getline(file, line)){
        if(has_start_codon(line)){
            gseq = line;
        }
        else{
            gseq = "ATG"+line;
        }
        if(!has_stop_codon(gseq)){
            gseq = gseq+"TAA";
        }
        hgt_pool.push_back(gseq);

    }
    file.close();
    std::cout<<"HGT pool has "<<hgt_pool.size()<<" genes\n";
    if(hgt_pool.size() > 1){
        std::cout<<hgt_pool[0].substr(0,50)<<"\n";
        std::cout<<hgt_pool[1].substr(0,50)<<"\n";
    }
    std::cout<<"----------------------------------------\n";
    std::cout<<"Shuffling HGT pool...\n";
    //std::srand ( unsigned ( std::time(0) ) );


    // std::srand ( unsigned ( config.rand_seed ) );
    
    std::random_shuffle(hgt_pool.begin(), hgt_pool.end());
    
    if(hgt_pool.size() > 1){
        std::cout<<hgt_pool[0].substr(0,50)<<"\n";
        std::cout<<hgt_pool[1].substr(0,50)<<"\n";
    }
    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading population tree...\n";
    std::cout<<config.itree<<"\n";
    std::map<int, std::vector<int> > tchilds;
    std::map<int,int> genome_parents;
    std::ifstream tfile(config.itree);
    int tn, tp;
    while( tfile >>  tn){
        tfile >> tp;
        //std::cout<<tn<<" "<<tp<<"\n";
        if( tchilds.find(tp) == tchilds.end() ){
            tchilds[tp] = std::vector<int>();
        }
        tchilds[tp].push_back(tn);
        genome_parents[tn] = tp;
    }
    tfile.close();
    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading substitution matrix...\n";
    std::cout<<config.isubm<<"\n";
    char **subm = read_subsitution_matrix(config.isubm);

    for(int i=0; i<4; i++){
        std::cout<<int2nuc[i]<<": ";
        for(int j=0; j<100; j++){
            std::cout<<subm[i][j];
        }
        std::cout<<"\n";
    }
    

    std::cout<<"----------------------------------------\n";
    int NOF_GENOMES = genome_parents.size();


    std::cout<<"Evolving "<<genome_parents.size()<<" genomes ...\n";


    //std::map<int,int> genome_parents;

    //std::map<GeneID,GeneID> gene_parents;
    std::map< std::pair<int,int>, std::pair<int,int> > gene_parents;

    // for chimeras records. this vector keeps track of all the chimeras that
    //  are generated during the evolution.
    //  the key of the map is a pair of integers (acceptor_genome_id,
    //  acceptor_gene_id) that identifies the gene to which the contribution
    // is added. the value is a ChimeraRecord struct that contains the list
    // of contributions that are added to the acceptor gene.
    // std::map<std::pair<int,int>, ChimeraRecord> chimeras;
    
    ChimeraLog chimera_log;

    //std::vector<Genome*> genomes;
    std::map<int, Genome*> genomes;

    //genomes.push_back(root_genome);
    genomes[0] = root_genome;
    //genome_parents[0] = -1;
    for(Locus &l : root_genome->loci){
        gene_parents[ std::pair<int,int>(0, l.id) ] = std::pair<int,int>(-1, -1);
    }

    int global_gene_id = root_genome->loci.size();

    //for(int run =0; run < NOF_GENOMES; run++){
    std::queue<int> treequeue;
    treequeue.push(0);
    if(!treequeue.empty()){
        int tn = treequeue.front(); treequeue.pop();

        std::map<int, std::vector<int> >::iterator tit = tchilds.find(tn);
        if(tit != tchilds.end()){
            for(int& ttn : tit->second){
                treequeue.push(ttn);
            }
        }

    }



    for(Locus &l : root_genome->loci){
        if(l.end >= root_genome->sequence.size()){
            std::cout<<"invalid loci end on root genome "<<l.end<<" "<<root_genome->sequence.size()<<"\n";
            exit(1);
        }
    }


    int globalcount = 0;
    while(!treequeue.empty()){
        int tn = treequeue.front(); treequeue.pop();
        std::cout<<"Generating genome "<<tn<<" ["<<globalcount<<"]\n";
        globalcount++;


        std::map<int, std::vector<int> >::iterator tit = tchilds.find(tn);
        if(tit != tchilds.end()){
            for(int& ttn : tit->second){
                treequeue.push(ttn);
            }
        }

        int current_genome_id = tn;

        //std::cout<<"genome "<<tn<<"\n";

        //int parent_genome_id = randint(0, genomes.size());
        int parent_genome_id = genome_parents[current_genome_id];

        Genome *parent_genome = genomes[parent_genome_id];
        std::cout<<"parent genome "<<parent_genome_id<<"\n";
        Genome *new_genome = parent_genome->clone();

        //genome_parents[run+1] = parent_genome_id;
        for(int i=0; i<new_genome->loci.size(); i++){
        //    gene_parents[ GeneID(run+1, &(new_genome->loci[i])) ] = GeneID(parent_genome_id, &(parent_genome->loci[i]));
            gene_parents[std::pair<int,int>(current_genome_id, new_genome->loci[i].id)] = std::pair<int,int>(parent_genome_id, parent_genome->loci[i].id);
        }

        for(Locus &l : new_genome->loci){
            if(l.end >= new_genome->sequence.size()){
                std::cout<<"invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                exit(1);
            }
        }

        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------


        for(int gv = 0; gv< parent_genome->loci.size(); gv++){
            if( generate_unif() <= config.gene_variation_prob){
                
#ifdef VERBOSE
                std::cout<<"--------------------\n";
                std::cout<<"altering gene "<<gv<<"\n";
#endif

                bool *constrained = (bool*)calloc(new_genome->sequence.size(), sizeof(bool));
                for(Locus &l : new_genome->loci){
                    if((l.start >= new_genome->loci[gv].start-3) && (l.start <= new_genome->loci[gv].end)){
                        for(int i=l.start; i<l.start+3; i++){
                            constrained[i]  =true;
                        }
                    }
                    if((l.end <= new_genome->loci[gv].start-3) && (l.end <= new_genome->loci[gv].end)){
                        for(int i=l.end-3; i<l.end; i++){
                            constrained[i]  =true;
                        }
                    }
                }

                int current_start = new_genome->loci[gv].start;
                int gene_length = new_genome->loci[gv].end - new_genome->loci[gv].start;
                int new_gene_length = gene_length;
                int current_shift = 0;
                int total_altered = 0;
                for(int p=0; p<gene_length; p++){
                    if(!constrained[p]){
                        
                        if( generate_unif() <= config.locus_variation_prob){
                            total_altered++;
                            int alteration = randint(0,4);
                            if(alteration == 0){
                                //snp
                                //std::cout<<new_genome->sequence.size()<<"["<<(current_start + p + current_shift)<<"]"<<current_start<<" "<<p<<" "<<current_shift<<"\n";
                                
                                //new_genome->sequence[current_start + p + current_shift] = int2nuc[randint(0,4)];

                                if(current_start + p + current_shift + 3 < new_genome->sequence.size()){
                                    //std::string cod = int2codons[randint(0, int2codons_len)];
                                    //new_genome->sequence[current_start + p + current_shift] = cod[0];
                                    //new_genome->sequence[current_start + p + current_shift +1] = cod[1];
                                    //new_genome->sequence[current_start + p + current_shift +2] = cod[2];

                                    do{
                                        //new_genome->sequence[current_start + p + current_shift] = 
                                        //    substitute( new_genome->sequence[current_start + p + current_shift],subm);
                                        new_genome->sequence[current_start + p + current_shift +1] = 
                                            substitute( new_genome->sequence[current_start + p + current_shift +1],subm);
                                        //new_genome->sequence[current_start + p + current_shift +2] = 
                                        //    substitute( new_genome->sequence[current_start + p + current_shift +2],subm);
                                    }while(to_revert(new_genome->sequence, current_start + p + current_shift ));
                                }
                            }
                            else if(alteration == 1){
                                //delete
                                if(current_start + p + current_shift + 3 < new_genome->sequence.size()){
                                    new_genome->sequence = 
                                        new_genome->sequence.substr(0, current_start + p + current_shift)
                                        +
                                        new_genome->sequence.substr(current_start + p + current_shift + 3)
                                        ;

                                    current_shift -= 3;

                                    for(Locus &l : new_genome->loci){
                                        if(l.start >= current_start +p){
                                            l.start -= 3;
                                        }
                                        if(l.end >= current_start +p){
                                            l.end -= 3;
                                        }
                                    }
                                }
                                // else{
                                //     new_genome->sequence = 
                                //         new_genome->sequence.substr(0, current_start + p + current_shift)
                                //         ;
                                // }
                                // current_shift -= 1;

                                // for(Locus &l : new_genome->loci){
                                //     if(l.start >= current_start +p){
                                //         l.start -= 1;
                                //     }
                                //     if(l.end >= current_start +p){
                                //         l.end -= 1;
                                //     }
                                // }
                            }
                            else{
                                //insert
                                if(current_start + p + current_shift + 3 < new_genome->sequence.size()){
                                    std::string cod = int2codons[randint(0, int2codons_len)];

                                    new_genome->sequence = 
                                        new_genome->sequence.substr(0, current_start + p + current_shift)
                                        +
                                        cod[0]
                                        +
                                        cod[1]
                                        +
                                        cod[2]
                                        +
                                        new_genome->sequence.substr(current_start + p + current_shift)
                                        ;

                                    current_shift += 3;

                                    for(Locus &l : new_genome->loci){
                                        if(l.start > current_start +p){
                                            l.start += 3;
                                        }
                                        if(l.end > current_start +p){
                                            l.end += 3;
                                        }
                                    }
                                }
                                // else{
                                //     new_genome->sequence = 
                                //         new_genome->sequence + int2nuc[randint(0,4)]
                                //         ;
                                // }
                                // current_shift += 1;

                                // for(Locus &l : new_genome->loci){
                                //     if(l.start > current_start +p){
                                //         l.start += 1;
                                //     }
                                //     if(l.end > current_start +p){
                                //         l.end += 1;
                                //     }
                                // }
                            }
                        }
                    }
                }
                new_gene_length += current_shift;
                free(constrained);

#ifdef VERBOSE
                std::cout<<"old gene length "<< gene_length <<"; altered "<<total_altered<<" positions; relative new length "<<current_shift<<"; new length "<<new_gene_length<<"\n";
                std::cout<<"result "<<new_genome->sequence.size()<< " "<<new_genome->loci.back()<<"\n";
#endif
                for(Locus &l : new_genome->loci){
                    if(l.end >= new_genome->sequence.size()){
                        std::cout<<"invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                        exit(1);
                    }
                }
                
            }
            if( generate_unif() <= config.gene_duplication_prob ){
#ifdef VERBOSE
                std::cout<<"duplicating "<<new_genome->loci[gv]<<"\n";
#endif

                std::string new_sequence = new_genome->sequence.substr(  new_genome->loci[gv].start,  new_genome->loci[gv].end - new_genome->loci[gv].start );
                std::pair<int, std::string*> new_gene = std::pair<int, std::string*>(-1, &new_sequence);

                //Locus *parent_gene_id = &(new_genome->loci[gv]), *new_gene_id;
                int parent_gene_id  = new_genome->loci[gv].id;
                int new_gene_id = global_gene_id;
                global_gene_id++;


                gene_parents[ std::pair<int,int>(current_genome_id, new_gene_id) ] = std::pair<int,int>(current_genome_id,parent_gene_id);
                    

                int strand = randint(2); if(strand==0) strand = -1;

                
                int start_position = randint( new_genome->sequence.size() );
                if((start_position == 0) || (start_position == new_genome->sequence.size()-1)){
                    if(start_position == 0){
                        new_genome->sequence = (*new_gene.second) + new_genome->sequence;
                        for(Locus &l  : new_genome->loci){
                            l.start +=  new_gene.second->size();
                            l.end +=  new_gene.second->size();
                        }
                        new_genome->loci.insert( new_genome->loci.begin(), Locus(new_gene_id, 0, new_gene.second->size(), strand) );
                        
                        //new_gene_id = &(new_genome->loci[0]);
                    }
                    else{
                        new_genome->loci.push_back(Locus(new_gene_id, new_genome->sequence.size(), new_genome->sequence.size()+new_gene.second->size(), strand));
                        new_genome->sequence =  new_genome->sequence + (*new_gene.second);

                        //new_gene_id = &(new_genome->loci[new_genome->loci.size()-1]);
                    }
                }
                else{
                    bool *covered = (bool*)calloc( new_genome->sequence.size(), sizeof(bool));
                    for(Locus &l : new_genome->loci){
                        for(int i=l.start; i<l.end; i++){
                            covered[i] = true;
                        }
                    }

                    while(start_position > 0){
                        if(!covered[start_position]){
                            break;
                        }
                        start_position--;
                    }
#ifdef VERBOSE
                    std::cout<<"adding at position "<<start_position<<"\n";
#endif

                    if(start_position == 0){
                        new_genome->sequence = (*new_gene.second) + new_genome->sequence;
                        for(Locus &l  : new_genome->loci){
                            l.start +=  new_gene.second->size();
                            l.end +=  new_gene.second->size();
                        }
                        new_genome->loci.insert( new_genome->loci.begin(), Locus(new_gene_id, 0, new_gene.second->size(), strand) );

                        //new_gene_id = &(new_genome->loci[0]);
                    }
                    else{
                        new_genome->sequence = 
                            new_genome->sequence.substr(0, start_position) + 
                            (*new_gene.second) 
                            //+ new_genome->sequence.substr(start_position, new_genome->sequence.size() - start_position);
                            + new_genome->sequence.substr(start_position);
                        for(Locus &l  : new_genome->loci){
                            if(l.start >= start_position){
                                l.start +=  new_gene.second->size();
                                l.end +=  new_gene.second->size();
                            }
                        }

                        new_genome->loci.push_back( Locus(new_gene_id, start_position, start_position + new_gene.second->size(), strand));
                        
                        //new_gene_id = &(new_genome->loci[new_genome->loci.size()-1]);

                        std::sort(new_genome->loci.begin(), new_genome->loci.end());
                    }
                    free(covered);
                }

                //gene_parents[ GeneID(run+1, new_gene_id) ] = GeneID(parent_genome_id,parent_gene_id);

                for(Locus &l : new_genome->loci){
                    if(l.end >= new_genome->sequence.size()){
                        std::cout<<"invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                        exit(1);
                    }
                }
             }

        }

        for(Locus &l : new_genome->loci){
            if(l.end >= new_genome->sequence.size()){
                std::cout<<"invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                exit(1);
            }
        }


        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        int deleted_genes = 0;
        int deleted_nucleotides = 0;
        std::vector< std::pair<int, std::string*> > genes_to_add;

#ifdef VERBOSE
        std::cout<<"gene variation set size is "<<std::ceil(  parent_genome->loci.size() * config.geneset_variation )<<" / "<<parent_genome->loci.size()<<"\n";
#endif
        for(int gv = 0; gv<std::ceil(  parent_genome->loci.size() * config.geneset_variation ); gv++){
            if( generate_unif() <= config.geneset_variation_remove){
                deleted_genes++;
                int gene_to_delete = randint( new_genome->loci.size() );

#ifdef VERBOSE
                std::cout<<"deleting "<<new_genome->loci[gene_to_delete]<<"\n";
#endif

                bool *covered = (bool*)calloc( new_genome->sequence.size(), sizeof(bool));
                for(int i=0; i<new_genome->loci.size(); i++){
                    if(i != gene_to_delete){
                        Locus l = new_genome->loci[i];
                        for(int i=l.start; i<l.end; i++){
                            covered[i] = true;
                        }
                    }
                }

                int gene_to_delete_start = new_genome->loci[gene_to_delete].start;
                int gene_to_delete_end = new_genome->loci[gene_to_delete].end;
                int gene_to_delete_id = new_genome->loci[gene_to_delete].id;

                int range_start = gene_to_delete_start;
                int range_length;
                int c_deleted_nuc = 0;

                for(int i=range_start+1; i<gene_to_delete_end; i++){
                    if(covered[i] != covered[i-1]){
                        if(covered[i]){
                            //we are closing a uncovering island
                            range_length = i - range_start;

                            for(Locus &l : new_genome->loci){
                                if(l.start + c_deleted_nuc >= i ){
                                    l.start -= range_length;
                                }
                                if(l.end + c_deleted_nuc >= i){
                                    l.end -= range_length;
                                }
                            }

                            new_genome->sequence = 
                                new_genome->sequence.substr(0, range_start - c_deleted_nuc) + 
                                new_genome->sequence.substr(i - c_deleted_nuc);

                            for(Locus &l : new_genome->loci){
                                if(l.end >= new_genome->sequence.size()){
                                    std::cout<<"(internal1) invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                                    exit(1);
                                }
                            }

                            range_start = -1;

                            deleted_nucleotides += range_length;
                            c_deleted_nuc += range_length;
                        }
                        else{
                            //we are opening a uncovering island
                            range_start = i;
                        }
                    }
                }
                if(range_start != -1){
                    int i = gene_to_delete_end;
                    range_length = i - range_start;

                    for(Locus &l : new_genome->loci){
                        if(l.start + c_deleted_nuc >= i){
                            l.start -= range_length;
                        }
                        if(l.end + c_deleted_nuc >= i ){
                            l.end -= range_length;
                        }
                    }

                    new_genome->sequence = 
                        new_genome->sequence.substr(0, range_start - c_deleted_nuc) + 
                        new_genome->sequence.substr(i - c_deleted_nuc);

                    for(Locus &l : new_genome->loci){
                        if(l.end >= new_genome->sequence.size()){
                            std::cout<<"(internal2) invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                            exit(1);
                        }
                    }

                    deleted_nucleotides += range_length;
                }


                free(covered);

                new_genome->loci.erase( new_genome->loci.begin() + gene_to_delete );

#ifdef VERBOSE
                std::cout<<"genome size is now "<<new_genome->sequence.size()<<"\n";
                std::cout<<"last locus is "<<new_genome->loci[ new_genome->loci.size()-1 ]<<"\n";
#endif

                for(Locus &l : new_genome->loci){
                    if(l.end >= new_genome->sequence.size()){
                        std::cout<<"invalid loci end  "<<l.end<<" "<<new_genome->sequence.size()<<"\n";
                        exit(1);
                    }
                }

            }
            if( generate_unif() <= config.geneset_variation_add){
                std::pair<int, std::string*> new_gene = generate_new_gene(hgt_pool, min_root_gene_length, max_root_gene_length);
                genes_to_add.push_back(new_gene);

                int strand = randint(2); if(strand==0) strand = -1;

                
                int parent_gene_id  = -1;
                int new_gene_id = global_gene_id;
                global_gene_id++;
                gene_parents[ std::pair<int,int>(current_genome_id, new_gene_id) ] = std::pair<int,int>(-1,parent_gene_id);

                
                int start_position = randint( new_genome->sequence.size() );
                if((start_position == 0) || (start_position == new_genome->sequence.size()-1)){
                    if(start_position == 0){
                        new_genome->sequence = (*new_gene.second) + new_genome->sequence;
                        for(Locus &l  : new_genome->loci){
                            l.start +=  new_gene.second->size();
                            l.end +=  new_gene.second->size();
                        }
                        new_genome->loci.insert( new_genome->loci.begin(), Locus(new_gene_id, 0, new_gene.second->size(), strand) );

                        //gene_parents[ GeneID(run+1, &(new_genome->loci[0])) ] = GeneID(-1,NULL);
#ifdef VERBOSE
                        std::cout<<"adding at position 0 a gene "<<new_genome->loci[0]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
#endif
                    }
                    else{
                        new_genome->loci.push_back(Locus(new_gene_id, new_genome->sequence.size(), new_genome->sequence.size()+new_gene.second->size(), strand));
                        new_genome->sequence =  new_genome->sequence + (*new_gene.second);
                        //gene_parents[ GeneID(run+1, &(new_genome->loci[ new_genome->loci.size()-1 ])) ] = GeneID(-1,NULL);
#ifdef VERBOSE
                        std::cout<<"adding at position last a gene "<<new_genome->loci[ new_genome->loci.size()-1 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
#endif
                    }


                }
                else{
                    std::cout<<"size "<< new_genome->sequence.size() <<"\n";
                    bool *covered = (bool*)calloc( new_genome->sequence.size(), sizeof(bool));
                    if(covered==NULL){
                        std::cout<<"OPSSSSSSSSSSSSSSS\n";
                    }
                    for(Locus &l : new_genome->loci){
                        for(int i=l.start; i<l.end; i++){
                            if(i>=new_genome->sequence.size() ){
                                std::cout<<"OPSSSSSSSSSSSSSSS i "<<i<<" "<< new_genome->sequence.size() <<"\n";
                            }
                            covered[i] = true;
                        }
                    }

                    while(start_position > 0){
                        if(!covered[start_position]){
                            break;
                        }
                        start_position--;
                    }
                    //std::cout<<"adding at position "<<start_position<<" a gene with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";

                    if(start_position == 0){
                        new_genome->sequence = (*new_gene.second) + new_genome->sequence;
                        for(Locus &l  : new_genome->loci){
                            l.start +=  new_gene.second->size();
                            l.end +=  new_gene.second->size();
                        }
                        new_genome->loci.insert( new_genome->loci.begin(), Locus(new_gene_id, 0, new_gene.second->size(), strand) );

                        //gene_parents[ GeneID(run+1, &(new_genome->loci[0])) ] = GeneID(-1,NULL);
#ifdef VERBOSE
                        std::cout<<"adding at position "<<start_position<<" a gene "<<new_genome->loci[ 0 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
#endif
                    }
                    else{
                        new_genome->sequence = 
                            new_genome->sequence.substr(0, start_position) + 
                            (*new_gene.second) 
                            //+ new_genome->sequence.substr(start_position, new_genome->sequence.size() - start_position);
                            + new_genome->sequence.substr(start_position);
                        for(Locus &l  : new_genome->loci){
                            if(l.start >= start_position){
                                l.start +=  new_gene.second->size();
                                l.end +=  new_gene.second->size();
                            }
                        }
                        new_genome->loci.push_back( Locus(new_gene_id, start_position, start_position + new_gene.second->size(), strand));
                        
                        //gene_parents[ GeneID(run+1, &(new_genome->loci[ new_genome->loci.size()-1 ])) ] = GeneID(-1,NULL);
#ifdef VERBOSE
                        std::cout<<"adding at position "<<start_position<<" a gene "<<new_genome->loci[ new_genome->loci.size()-1 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
#endif

                        std::sort(new_genome->loci.begin(), new_genome->loci.end());
                    }
                    free(covered);
                }

                if(new_gene.first == -1){
                    delete new_gene.second;
                }
                else{
                    hgt_pool.erase( hgt_pool.begin() + new_gene.first );
                }
            }
        }
        
        std::cout<<deleted_genes<<" genes deleted\n";
        std::cout<<deleted_nucleotides<<" deleted nucleotides\n";

        std::cout<<genes_to_add.size()<<" genes added\n";

        std::cout<<"the new genome has a total of "<<new_genome->loci.size()<<" genetic loci\n";
        std::cout<<"the length of the new genome is "<<new_genome->sequence.size()<<" nucleotides\n";


        // ---------------------
        // ---------------------
        // Gene fusion section beginning
        // ---------------------
        // ---------------------

        if( generate_unif() <= config.gene_fusion_prob ){
#ifdef VERBOSE
std::cout << "fusing in genome " << current_genome_id << "\n";
#endif
            // ---------------------
            // duplication - max 1 duplication event per genome
            // ---------------------
            if (generate_unif() <= config.sub_gene_duplication_prob) {

                if (new_genome->loci.size() == 0) {
                    // no genes to duplicate, skipping
                    std::cout<<"no genes to duplicate, skipping\n";
                } else {
                    // random source gene
                    int source_gene_index = randint(new_genome->loci.size());

                    const Locus source_gene = new_genome->loci[source_gene_index];
                    int source_gene_len = source_gene.end - source_gene.start;
                    if (source_gene_len < 3 ) {
                        // gene to short
                        std::cout<<"source gene too short for sub-gene duplication, skipping\n";
                    } else {

                        //  i must keep x3 pattern
                        int source_gene_total_codons = source_gene_len / 3;
                        int source_gene_min_codon_start = 1;
                        int source_gene_max_codon_start = source_gene_total_codons - 1;

                        if (source_gene_max_codon_start - source_gene_min_codon_start < 1) {
                            // gene too short for sub-gene duplication, skipping
                            std::cout<<"source gene too short for sub-gene duplication, skipping\n";
                        
                        } else {
                        
                            int source_gene_codon1 = randint(source_gene_min_codon_start, source_gene_max_codon_start); // [first..last)
                            int source_gene_codon2 = randint(source_gene_codon1 + 1, source_gene_max_codon_start + 1); // [c1+1..last]

                            int source_gene_range_offset = 3 * source_gene_codon1;
                            int source_gene_range_length = 3 * (source_gene_codon2 - source_gene_codon1);
                            int source_gene_absolute_start = source_gene.start + source_gene_range_offset;

                            std::string source_gene_sub_seq = new_genome->sequence.substr(source_gene_absolute_start, source_gene_range_length);


                            // dice if same gene or not
                            if (generate_unif() <= config.intrasub_gene_duplication_prob) { // same gene

                                // concat to the end of duplication site
#ifdef VERBOSE
std::cout << "duplicating a sequence of length " << source_gene_range_length << " from gene " << source_gene.id << " at offset " << source_gene_range_offset << " to the same gene\n";
#endif
                                new_genome->sequence = 
                                    new_genome->sequence.substr(0, source_gene_absolute_start + source_gene_range_length) +
                                    source_gene_sub_seq +
                                    new_genome->sequence.substr(source_gene_absolute_start + source_gene_range_length);


                                // update loci positions
                                for (Locus& l : new_genome->loci) {
                                    if (l.end <= source_gene_absolute_start + source_gene_range_length) {
                                        // skip: gene before the insertion point, no shift
                                    } else if (l.start >= source_gene_absolute_start + source_gene_range_length) {
                                        // gene after the insertion point: shift both start and end
                                        l.start += source_gene_range_length;
                                        l.end   += source_gene_range_length;
                                    } else {
                                        // gene that crosses the insertion point (= the source): only end advances
                                        l.end += source_gene_range_length;
                                    }
                                }

                                
                                ChimeraContribution chimera_contrib = make_chimera_contribution(
                                    std::move(make_locus(
                                        current_genome_id,
                                        source_gene.id
                                    )),
                                    source_gene_range_offset,
                                    source_gene_range_length
                                    #ifdef DEBUG
                                    , source_gene_sub_seq
                                    #endif
                                );
#ifdef VERBOSE
std::cout << "recording chimera contribution: source gene " << source_gene.id << " (offset " << source_gene_range_offset << ", length " << source_gene_range_length << ") contributes to the same gene " << source_gene.id << " at offset " << source_gene_absolute_start + source_gene_range_length << "\n";
std::cout << "dump of the chimera contribution:\n" << chimera_contrib << "\n";
#endif
                                std::vector<ChimeraContribution> contributions = {chimera_contrib};
                                
                                ChimeraRecord chimera_record = make_chimera_record(
                                    ChimeraEventType::GENE_FUSION_INTRA_SUB_GENE_DUPLICATION,
                                    std::move(
                                        make_chimera_acceptor(
                                            std::move(
                                                make_locus(
                                                current_genome_id,
                                                source_gene.id
                                                )
                                            ),
                                            source_gene_absolute_start + source_gene_range_length
                                        )
                                    ),
                                    std::move(contributions)
                                );

#ifdef VERBOSE
std::cout << "dump of the chimera record:\n" << chimera_record << "\n";
#endif
                                
                                chimera_log.add_chimera_event(std::move(chimera_record));
#ifdef DEBUG
chimera_contrib.contribution_sequence = source_gene_sub_seq;
#endif

                            } else if (generate_unif() <= config.intersub_gene_duplication_prob) { // different gene

                                // select a random gene and a random position in it, and insert the duplicated sequence there

                                // generate a random gene target index, if it is the same gene as the source,
                                // the index is re-generated until a different gene is selected.
                                int target_gene_index = source_gene_index;
                                while (target_gene_index == source_gene_index) {
                                    target_gene_index = randint(new_genome->loci.size());
                                    if (target_gene_index == source_gene_index) {
                                        std::cout << "same gene selected for duplication, re-rolling\n";
                                    }
                                }

                                // select a random position in the target gene and insert the duplicated sequence there
                                Locus& target_gene = new_genome->loci[target_gene_index];
                                int target_gene_len = target_gene.end - target_gene.start;
                                int target_total_codons = target_gene_len / 3;
                                if (target_total_codons < 2) {
                                    // target to short
                                    std::cout << "target gene too short, skipping\n";
                                } else {
                                    int target_gene_insert_position = randint(1, target_total_codons); // [1..target_total_codons-1]
                                    int target_gene_insert_offset = 3 * target_gene_insert_position;
                                    int target_gene_insert_pos = target_gene.start + target_gene_insert_offset;

                                    new_genome->sequence =
                                        new_genome->sequence.substr(0, target_gene_insert_pos) +
                                        source_gene_sub_seq +
                                        new_genome->sequence.substr(target_gene_insert_pos);

                                    // update loci positions
                                    for (Locus& l : new_genome->loci) {
                                        if (l.end <= target_gene_insert_pos) {
                                            // skip: gene before the insertion point, no shift
                                        } else if (l.start >= target_gene_insert_pos) {
                                            // gene after the insertion point: shift both start and end
                                            l.start += source_gene_range_length;
                                            l.end   += source_gene_range_length;
                                        } else {
                                            // gene that crosses the insertion point (= the target): only end advances
                                            l.end += source_gene_range_length;
                                        }
                                    }

                                    // tracking contribution to chimera
                                    ChimeraContribution chimera_contrib = make_chimera_contribution(
                                        std::move(make_locus(
                                            current_genome_id,
                                            source_gene.id
                                        )),
                                        source_gene_range_offset,
                                        source_gene_range_length
                                        #ifdef DEBUG
                                        , source_gene_sub_seq
                                        #endif
                                    );
#ifdef VERBOSE
std::cout << "recording chimera contribution: source gene " << source_gene.id << " (offset " << source_gene_range_offset << ", length " << source_gene_range_length << ") contributes to target gene " << target_gene.id << " at offset " << target_gene_insert_pos << "\n";
std::cout << "dump of the chimera contribution:\n" << chimera_contrib << "\n";
#endif
                                    std::vector<ChimeraContribution> contributions = {chimera_contrib};
                                    
                                    ChimeraRecord chimera_record = make_chimera_record(
                                        ChimeraEventType::GENE_FUSION_INTER_SUB_GENE_DUPLICATION,
                                        std::move(
                                            make_chimera_acceptor(
                                                std::move(
                                                    make_locus(
                                                        current_genome_id,
                                                        target_gene.id
                                                    )
                                                ),
                                                target_gene_insert_pos
                                            )
                                        ),
                                        std::move(contributions)
                                    );

#ifdef VERBOSE
std::cout << "dump of the chimera record:\n" << chimera_record << "\n";
#endif
                                    chimera_log.add_chimera_event(std::move(chimera_record));

                                    #ifdef VERBOSE
                                    std::cout << "recording chimera contribution: source gene " << source_gene.id << " (offset " << source_gene_range_offset << ", length " << source_gene_range_length <<
                                        ") contributes to target gene " << target_gene.id << " at offset " << target_gene_insert_pos << "\n";
                                    #endif
                                    #ifdef DEBUG
                                    chimera_contrib.contribution_sequence = source_gene_sub_seq;
                                    #endif

                                    // check if the new genome is valid
                                    for (Locus& l : new_genome->loci) {
                                        if (l.end >= (int)new_genome->sequence.size()) {
                                            std::cout << "(chimera dup) invalid loci end "
                                                        << l.end << " " << new_genome->sequence.size() << "\n";
                                            exit(1);
                                        }
                                    }

                                }
                            } else { // no duplication, skipping
                                std::cout << "no duplication: duplicated sequence dropped, skipping\n";
                            }

                        }
                    }
                }
            }


            // ---------------------
            // extended deletion - max 1 extended deletion event per genome
            // ---------------------

            // TODO: other fusion events
            // extended deletion
            // if (generate_unif() <= config.extended_gene_deletion_prob) {
                // TODO add keep prob -> if, drop otherwise
            // }

        }

        // ---------------------
        // ---------------------
        // Gene fission section beginning
        // ---------------------
        // ---------------------
        // TODO

        // add new genome to the population
        genomes[current_genome_id] = new_genome;

        std::cout<<"----------------------------------------\n";
    }

    std::cout<<"HGT pool has "<<hgt_pool.size()<<" genes\n";

    std::cout<<"----------------------------------------\n";

    std::string oprefix(config.output_prefix);

    //std::map<int,int> genome_parents;
    //std::map<GeneID,GeneID> gene_parents;
    //std::vector<Genome*> genomes;

    std::ofstream myfile;

    std::cout<<"writing genomic sequences "<<oprefix<<".genome_sequences ...\n";
    myfile.open (oprefix + ".genome_sequences");
    int i=0;
    //for(auto genome : genomes){
    for(auto const& x : genomes){
        Genome *genome = x.second;
        myfile<<i<<" "<<genome->sequence<<"\n";
        i++;
    }
    myfile.flush();
    myfile.close();

    /*std::cout<<"writing genomic parenthood "<<oprefix<<".genome_parents ...\n";
    //std::ofstream myfile;
    myfile.open (oprefix + ".genome_parents");
    for(auto const &x : genome_parents){
        myfile<<x.first<<" "<<x.second<<"\n";
    }
    myfile.flush();
    myfile.close();*/

    std::cout<<"writing genetic parenthood "<<oprefix<<".gene_parents ...\n";
    //std::ofstream myfile;
    myfile.open (oprefix + ".gene_parents");
    for(auto const &x : gene_parents){
        //myfile<<x.first.first<<" "<<x.first.second<<" "<<x.second.first<<" "<<x.second.second<<"\n";
        myfile<<x.first.first<<" "<<x.first.second<<" "<<x.second.first<<" "<<x.second.second<<"\n";
    }
    myfile.flush();
    myfile.close();


    std::cout<<"writing genetic information "<<oprefix<<".genes ...\n";
    //std::ofstream myfile;
    myfile.open (oprefix + ".genes");
    //for(int gid=0; gid<genomes.size(); gid++){
    for(auto const& x :genomes){
        int gid = x.first;
        for(Locus &l : genomes[gid]->loci){
            if(l.strand == 1){
                myfile<<gid<<":"<<l.id<<":("<<l.start<<","<<(l.end-1)<<","<<l.strand<<") "<< genomes[gid]->sequence.substr(l.start, l.end-l.start)  <<"\n";
            }
            else{

                std::string t = genomes[gid]->sequence.substr(l.start, l.end-l.start);
                std::string s = genomes[gid]->sequence.substr(l.start, l.end-l.start);
                for(int i=0; i<s.size(); i++){
                    s[i] = Genome::rc_symbol( t[ s.size()-1-i ] );
                }
                // std::string s = genomes[gid]->sequence.substr(l.start, l.end-l.start);
                // for(int i=0; i<s.size(); i++){
                //     s[i] = Genome::rc_symbol( genomes[gid]->sequence[ l.end-i] );
                // }
                myfile<<gid<<":"<<l.id<<":("<<l.start<<","<<(l.end-1)<<","<<l.strand<<") "<< s <<"\n";
            }
        }
    }
    myfile.flush();
    myfile.close();

    std::cout<<"writing chimera events "<<oprefix<<".chimeras.csv\n";
    chimera_log.write_to_csv(oprefix + ".chimeras.csv");
    std::cout<<"writing chimera events "<<oprefix<<".chimeras.tsv\n";
    chimera_log.write_to_tsv(oprefix + ".chimeras.tsv");
};
