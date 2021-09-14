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


#define NOF_GENOMES 999

#define GENE_VARIATION_PROB 0.5	//probability of variation when ancestor gene is aquired
#define LOCUS_VARIATION_PROB 0.01 // 0.05  // probability of variating a nucleotide in a variated sequence


#define GENE_DUPLICATION_PROB 0.001	//probability of duplicating a gene
#define GENESET_VARIATION 0.01	//percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
#define GENESET_VARIATION_ADD 0.9 //probability that the variation is a a gene add
#define GENESET_VARIATION_REMOVE 0.1 //probability hat the variation if a gene removal

#define RAND_SEED 123456789


void usage(std::string cmd){
    std::cout<<"Usage: "<<cmd<<" root_genome.peg hgt_pool.hgt\n";
}




class Locus{
public:
    int id;
    int start;
    int end;
    int strand;// 1, -1

    Locus(){
        this->id = -1;
        this->start = -1;
        this->end = -1;
        this->strand = 0;
    }

    /*Locus(int _start, int _end, int _strand)
        : start(_start), end(_end), strand(_strand)
    {
        this->id = -1;
    }*/

    Locus(int _id, int _start, int _end, int _strand)
        : id(_id), start(_start), end(_end), strand(_strand)
    {
    }

    Locus(const Locus &l){
        this->id = l.id;
        this->start = l.start;
        this->end = l.end;
        this->strand= l.strand;
    }

    bool operator < (const Locus& a) const
    {
        if(start == a.start){
            if(end == a.end){
                return strand > a.strand;
            }
            return end < a.end;
        }
        return start < a.start;
    }

    friend std::ostream& operator<<(std::ostream& os, Locus const & tc) {
        return os << "("<<tc.id<<","<<tc.start<<","<<tc.end<<","<<tc.strand<<")";
    }
};


class Genome{
public:
    std::string sequence;
    std::vector<Locus> loci;

    Genome(){
        sequence = "";
    }

    Genome(const Genome& o){
        this->sequence = o.sequence;
        for(Locus x : o.loci){
            this->loci.push_back(x);
        }
    }

    Genome*
    clone(){
        Genome *g = new Genome();
        g->sequence = this->sequence;
        for(Locus x : this->loci){
            g->loci.push_back(x);
        }
        return g;
    }

    static 
    Genome* read_from_file(std::string ifile) {
        ifile.erase(std::remove_if(ifile.begin(), ifile.end(), [](unsigned char x){return std::isspace(x);}), ifile.end());
        Genome *g = new Genome();

        std::ifstream file(ifile);
        std::string line; 
        
        if(std::getline(file, line)){
            g->sequence = std::string(line);
            std::transform(g->sequence.begin(), g->sequence.end(),g->sequence.begin(), ::toupper);
            //std::cout<<"[sequence end]\n";
        }
        else{
            //std::cout<<"[no sequence]\n";
            return NULL;
        }

        int max_locus_id = 0;
        int c;
        int state = 0;
        int start, end;
        char strand;
        while(file >> c){
            if(state == 0){
                start = c;
                state = 1;
            }
            else if(state == 1){
                end = c;
                state = 2;
            }
            else{
                //std::cout<<start<<" "<<end<<" "<<c<<"\n";
                g->loci.push_back( Locus(max_locus_id, start,end,c) );
                max_locus_id++;

                state = 0;
            }
        }
        //std::cout<<"[loci end]\n";
        file.close();

        std::sort(g->loci.begin(), g->loci.end());
        
        return g;
    };

    static
    char rc_symbol(char c){
        if(c == 'A') return 'T';
        if(c == 'T') return 'A';
        if(c == 'C') return 'G';
        if(c == 'G') return 'C';
        return 'N';
    }

    std::vector<std::string>*
    get_gene_sequences(){
        std::vector<std::string> *v = new std::vector<std::string>();
        for(auto & locus : this->loci){
            if(locus.strand == 1){
                std::string s = this->sequence.substr(locus.start, locus.end-locus.start+1);
                v->push_back(s);
            }
            else{
                std::string s = this->sequence.substr(locus.start, locus.end-locus.start+1);
                for(int i=0; i<s.size(); i++){
                    s[i] = Genome::rc_symbol( this->sequence[ locus.end+1-i] );
                }
                v->push_back(s);
            }
        }
        return v;
    };


};


/*
class GeneID{
public:
    int genomeID;
    Locus *gene_locus;

    GeneID(){
        this->genomeID = -1;
        this->gene_locus = NULL;
    }

    GeneID(int _gid, Locus *_locus){
        this->genomeID = _gid;
        this->gene_locus = _locus;
    }

    GeneID(const GeneID &g){
        this->genomeID = g.genomeID;
        this->gene_locus = g.gene_locus;
    }

    bool operator < (const GeneID& a) const
    {
        if(genomeID == a.genomeID){
            if(gene_locus == NULL && a.gene_locus == NULL){
                return true;
            }
            else if(gene_locus != NULL && a.gene_locus == NULL){
                return true;
            }
            else if(gene_locus == NULL && a.gene_locus == NULL){
                return false;
            }
            else{
                return gene_locus < a.gene_locus;
            }
        }
        return genomeID < a.genomeID;
    }

};
*/



//AUG GUG and UUG -> ATG GTG TTG
bool
has_start_codon(std::string& s){
    if(s.size() < 3) return false;
    if((s[0]=='A')&&(s[1]=='T')&&(s[2]=='G')) return true;
    if((s[0]=='G')&&(s[1]=='T')&&(s[2]=='G')) return true;
    if((s[0]=='T')&&(s[1]=='T')&&(s[2]=='G')) return true;
    return false;
}

//TAA, TGA, TAG
bool
has_stop_codon(std::string& s){
    if(s.size() < 3) return false;
    if((s[s.size()-3]=='T')&&(s[s.size()-2]=='A')&&(s[s.size()-1]=='A')) return true;
    if((s[s.size()-3]=='T')&&(s[s.size()-2]=='G')&&(s[s.size()-1]=='A')) return true;
    if((s[s.size()-3]=='T')&&(s[s.size()-2]=='A')&&(s[s.size()-1]=='G')) return true;
    return false;
}


int randint (int i) { return std::rand()%i;}
int randint (int i,int j) { return (std::rand()%(j-i))+i;}


static const char int2nuc[] = {'A','C','G','T'};

std::string*
make_random_sequence(int min_length, int max_length){

    int length = randint(min_length, max_length) + 6;
    char * c = new char[length + 1];
    c[0] = 'A'; c[1] = 'T'; c[2] = 'G';
    c[length-3] = 'T'; c[length-2] = 'A'; c[length-1] = 'A'; 
    c[length] = 0;

    for(int i=0; i<length; i++){
        c[i] = int2nuc[ randint(0,4) ];
    }

    std::string *s = new std::string(c);
    return s;
};


std::pair<int, std::string*>
generate_new_gene(std::vector<std::string> &hgt_pool, int min_length, int max_length){
    if(hgt_pool.size() == 0){
        return std::pair<int, std::string*>(-1, make_random_sequence(min_length,max_length));
    }
    else{
        int i = randint(hgt_pool.size());
        return std::pair<int, std::string*>(i, &(hgt_pool[i]));
    }
};


int main(int argc, char** argv){
    std::cout<<argc<<"\n";
    if(argc!=4){
        usage(argv[0]);
        return 0;
    }


    /*
    for(int i=0; i<10; i++){
        std::string *s = make_random_sequence(10,20);
        std::cout<<(*s)<<"\n";
        delete s;
    }
    */

    /*
    Genome a;
    a.sequence = "ATTG";
    a.loci.push_back(Locus(2,10,1));
    Genome b = a;
    b.sequence[0]='T';
    b.loci[0].start = 3;
    std::cout<<a.sequence<<"\n";
    std::cout<<b.sequence<<"\n";
    std::cout<<a.loci[0]<<"\n";
    std::cout<<b.loci[0]<<"\n";
    */

   /*
   std::vector<Locus> t_loci;
   t_loci.push_back( Locus(0,10,1) );
   t_loci.push_back( Locus(11,20,1) );

   for(Locus &l : t_loci){
       l.start += 10;
       l.end += 10;
   }
   for(Locus &l : t_loci){
       std::cout<<l<<"\n";
   }
   */


    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading root genome...\n";
    Genome *root_genome = Genome::read_from_file(argv[1]);

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
    }
    std::cout<<"gene lenghts are between "<<min_root_gene_length<<" and "<<max_root_gene_length<<" \n";

    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading HGT pool...\n";
    std::vector<std::string> hgt_pool;
    std::ifstream file(argv[2]);
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
    std::srand ( unsigned ( RAND_SEED ) );
    std::random_shuffle(hgt_pool.begin(), hgt_pool.end());
    if(hgt_pool.size() > 1){
        std::cout<<hgt_pool[0].substr(0,50)<<"\n";
        std::cout<<hgt_pool[1].substr(0,50)<<"\n";
    }
    std::cout<<"----------------------------------------\n";
    std::cout<<"Evolving "<<NOF_GENOMES<<" genomes ...\n";


    std::map<int,int> genome_parents;
    //std::map<GeneID,GeneID> gene_parents;
    std::map< std::pair<int,int>, std::pair<int,int> > gene_parents;

    std::vector<Genome*> genomes;

    genomes.push_back(root_genome);
    genome_parents[0] = -1;
    for(Locus &l : root_genome->loci){
        gene_parents[ std::pair<int,int>(0, l.id) ] = std::pair<int,int>(-1, -1);
    }


    std::mt19937_64 rng;
    rng.seed(RAND_SEED);
    std::uniform_real_distribution<double> unif(0, 1);

    int global_gene_id = root_genome->loci.size();


    for(int run =0; run < NOF_GENOMES; run++){
        std::cout<<"genome "<<run<<"\n";

        int parent_genome_id = randint(0, genomes.size());
        Genome *parent_genome = genomes[parent_genome_id];
        std::cout<<"parent genome "<<parent_genome_id<<"\n";
        Genome *new_genome = parent_genome->clone();

        genome_parents[run+1] = parent_genome_id;
        for(int i=0; i<new_genome->loci.size(); i++){
        //    gene_parents[ GeneID(run+1, &(new_genome->loci[i])) ] = GeneID(parent_genome_id, &(parent_genome->loci[i]));
            gene_parents[std::pair<int,int>(run+1, new_genome->loci[i].id)] = std::pair<int,int>(parent_genome_id, parent_genome->loci[i].id);
        }

        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------


        for(int gv = 0; gv< parent_genome->loci.size(); gv++){
            if( unif(rng) <= GENE_VARIATION_PROB){
                
                std::cout<<"--------------------\n";
                std::cout<<"altering gene "<<gv<<"\n";
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
                        if( unif(rng) <= LOCUS_VARIATION_PROB){
                            total_altered++;
                            int alteration = randint(0,4);
                            if(alteration == 0){
                                //snp
                                //std::cout<<new_genome->sequence.size()<<"["<<(current_start + p + current_shift)<<"]"<<current_start<<" "<<p<<" "<<current_shift<<"\n";
                                new_genome->sequence[current_start + p + current_shift] = int2nuc[randint(0,4)];
                            }
                            else if(alteration == 1){
                                //delete
                                if(current_start + p + current_shift + 1 < new_genome->sequence.size()){
                                    new_genome->sequence = 
                                        new_genome->sequence.substr(0, current_start + p + current_shift)
                                        +
                                        new_genome->sequence.substr(current_start + p + current_shift + 1)
                                        ;
                                }
                                else{
                                    new_genome->sequence = 
                                        new_genome->sequence.substr(0, current_start + p + current_shift)
                                        ;
                                }
                                current_shift -= 1;

                                for(Locus &l : new_genome->loci){
                                    if(l.start >= current_start +p){
                                        l.start -= 1;
                                    }
                                    if(l.end >= current_start +p){
                                        l.end -= 1;
                                    }
                                }
                            }
                            else{
                                //insert
                                if(current_start + p + current_shift + 1 < new_genome->sequence.size()){
                                    new_genome->sequence = 
                                        new_genome->sequence.substr(0, current_start + p + current_shift)
                                        +
                                        int2nuc[randint(0,4)]
                                        +
                                        new_genome->sequence.substr(current_start + p + current_shift)
                                        ;
                                }
                                else{
                                    new_genome->sequence = 
                                        new_genome->sequence + int2nuc[randint(0,4)]
                                        ;
                                }
                                current_shift += 1;

                                for(Locus &l : new_genome->loci){
                                    if(l.start > current_start +p){
                                        l.start += 1;
                                    }
                                    if(l.end > current_start +p){
                                        l.end += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                new_gene_length += current_shift;
                free(constrained);

                std::cout<<"old gene length "<< gene_length <<"; altered "<<total_altered<<" positions; relative new length "<<current_shift<<"; new length "<<new_gene_length<<"\n";
                std::cout<<"result "<<new_genome->sequence.size()<< " "<<new_genome->loci.back()<<"\n";
                
            }
            if( unif(rng) <= GENE_DUPLICATION_PROB ){
                std::cout<<"duplicating "<<new_genome->loci[gv]<<"\n";

                std::string new_sequence = new_genome->sequence.substr(  new_genome->loci[gv].start,  new_genome->loci[gv].end - new_genome->loci[gv].start );
                std::pair<int, std::string*> new_gene = std::pair<int, std::string*>(-1, &new_sequence);

                //Locus *parent_gene_id = &(new_genome->loci[gv]), *new_gene_id;
                int parent_gene_id  = new_genome->loci[gv].id;
                int new_gene_id = global_gene_id;
                global_gene_id++;


                gene_parents[ std::pair<int,int>(run+1, new_gene_id) ] = std::pair<int,int>(run+1,parent_gene_id);
                    

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
                    std::cout<<"adding at position "<<start_position<<"\n";

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
            }
        }


        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        int deleted_genes = 0;
        int deleted_nucleotides = 0;
        std::vector< std::pair<int, std::string*> > genes_to_add;

        std::cout<<"gene variation set size is "<<std::ceil(  parent_genome->loci.size() * GENESET_VARIATION )<<" / "<<parent_genome->loci.size()<<"\n";
        for(int gv = 0; gv<std::ceil(  parent_genome->loci.size() * GENESET_VARIATION ); gv++){
            if( unif(rng) <= GENESET_VARIATION_REMOVE){
                deleted_genes++;
                int gene_to_delete = randint( new_genome->loci.size() );

                std::cout<<"deleting "<<new_genome->loci[gene_to_delete]<<"\n";

                bool *covered = (bool*)calloc( new_genome->sequence.size(), sizeof(bool));
                for(int i=0; i<new_genome->loci.size(); i++){
                    if(i != gene_to_delete){
                        Locus l = new_genome->loci[i];
                        for(int i=l.start; i<l.end; i++){
                            covered[i] = true;
                        }
                    }
                }

                //int range_start = new_genome->loci[gene_to_delete].start;
                int range_start = new_genome->loci[gene_to_delete].start;
                int range_length;
                for(int i=range_start+1; i<new_genome->loci[gene_to_delete].end; i++){
                    if(covered[i] != covered[i-1]){
                        if(covered[i]){
                            //we are closing a covering island
                            range_length = i - range_start;

                            for(Locus &l : new_genome->loci){
                                if(l.start >= i){
                                    l.start -= range_length;
                                    l.end -= range_length;
                                }
                            }

                            new_genome->sequence = 
                                new_genome->sequence.substr(0, range_start) + 
                                new_genome->sequence.substr(i);

                            range_start = -1;

                            deleted_nucleotides += range_length;
                        }
                        else{
                            //we are opening a covering island
                            range_start = i;
                        }
                    }
                }
                if(range_start != -1){
                    int i = new_genome->loci[gene_to_delete].end;
                    range_length = i - range_start;

                    for(Locus &l : new_genome->loci){
                        if(l.start >= i){
                            l.start -= range_length;
                            l.end -= range_length;
                        }
                    }

                    new_genome->sequence = 
                        new_genome->sequence.substr(0, range_start) + 
                        new_genome->sequence.substr(i);

                    deleted_nucleotides += range_length;
                }


                free(covered);

                new_genome->loci.erase( new_genome->loci.begin() + gene_to_delete );

                std::cout<<"genome size is now "<<new_genome->sequence.size()<<"\n";
                std::cout<<"last locus is "<<new_genome->loci[ new_genome->loci.size()-1 ]<<"\n";

            }
            if( unif(rng) <= GENESET_VARIATION_ADD){
                std::pair<int, std::string*> new_gene = generate_new_gene(hgt_pool, min_root_gene_length, max_root_gene_length);
                genes_to_add.push_back(new_gene);

                int strand = randint(2); if(strand==0) strand = -1;

                
                int parent_gene_id  = -1;
                int new_gene_id = global_gene_id;
                global_gene_id++;
                gene_parents[ std::pair<int,int>(run+1, new_gene_id) ] = std::pair<int,int>(-1,parent_gene_id);

                
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
                        std::cout<<"adding at position 0 a gene "<<new_genome->loci[0]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
                    }
                    else{
                        new_genome->loci.push_back(Locus(new_gene_id, new_genome->sequence.size(), new_genome->sequence.size()+new_gene.second->size(), strand));
                        new_genome->sequence =  new_genome->sequence + (*new_gene.second);
                        //gene_parents[ GeneID(run+1, &(new_genome->loci[ new_genome->loci.size()-1 ])) ] = GeneID(-1,NULL);
                        std::cout<<"adding at position last a gene "<<new_genome->loci[ new_genome->loci.size()-1 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
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
                    //std::cout<<"adding at position "<<start_position<<" a gene with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";

                    if(start_position == 0){
                        new_genome->sequence = (*new_gene.second) + new_genome->sequence;
                        for(Locus &l  : new_genome->loci){
                            l.start +=  new_gene.second->size();
                            l.end +=  new_gene.second->size();
                        }
                        new_genome->loci.insert( new_genome->loci.begin(), Locus(new_gene_id, 0, new_gene.second->size(), strand) );

                        //gene_parents[ GeneID(run+1, &(new_genome->loci[0])) ] = GeneID(-1,NULL);
                        std::cout<<"adding at position "<<start_position<<" a gene "<<new_genome->loci[ 0 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";
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
                        std::cout<<"adding at position "<<start_position<<" a gene "<<new_genome->loci[ new_genome->loci.size()-1 ]<<" with "<<new_gene.second->size()<<" nucleotides; strand "<<strand<<"; hgt pool id "<<new_gene.first<<"\n";

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


                //if(new_gene.first != -1){
                //    hgt_pool.erase( hgt_pool.begin()+new_gene.first );
                //}

            }
        }
        
        std::cout<<deleted_genes<<" genes to be deleted\n";
        std::cout<<deleted_nucleotides<<" deleted nucleotides\n";

        std::cout<<genes_to_add.size()<<" genes added\n";

        std::cout<<"the new genome has a total of "<<new_genome->loci.size()<<" genetic loci\n";
        std::cout<<"the length of the new genome is "<<new_genome->sequence.size()<<" nucleotides\n";

        genomes.push_back(new_genome);

        std::cout<<"----------------------------------------\n";
    }

    std::cout<<"HGT pool has "<<hgt_pool.size()<<" genes\n";

    std::cout<<"----------------------------------------\n";

    std::string oprefix(argv[3]);

    //std::map<int,int> genome_parents;
    //std::map<GeneID,GeneID> gene_parents;
    //std::vector<Genome*> genomes;

    std::ofstream myfile;

    std::cout<<"writing genomic sequences "<<oprefix<<".genome_sequences ...\n";
    myfile.open (oprefix + ".genome_sequences");
    int i=0;
    for(auto genome : genomes){
        myfile<<i<<" "<<genome->sequence<<"\n";
        i++;
    }
    myfile.flush();
    myfile.close();

    std::cout<<"writing genomic parenthood "<<oprefix<<".genome_parents ...\n";
    //std::ofstream myfile;
    myfile.open (oprefix + ".genome_parents");
    for(auto const &x : genome_parents){
        myfile<<x.first<<" "<<x.second<<"\n";
    }
    myfile.flush();
    myfile.close();

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
    for(int gid=0; gid<genomes.size(); gid++){
        for(Locus &l : genomes[gid]->loci){
            if(l.strand == 1){
                myfile<<gid<<":"<<l.id<<":("<<l.start<<","<<l.end<<","<<l.strand<<") "<< genomes[gid]->sequence.substr(l.start, l.end-l.start)  <<"\n";
            }
            else{
                std::string s = genomes[gid]->sequence.substr(l.start, l.end-l.start+1);
                for(int i=0; i<s.size(); i++){
                    s[i] = Genome::rc_symbol( genomes[gid]->sequence[ l.end+1-i] );
                }
                myfile<<gid<<":"<<l.id<<":("<<l.start<<","<<l.end<<","<<l.strand<<") "<< s <<"\n";
            }

        }
    }
    myfile.flush();
    myfile.close();

};