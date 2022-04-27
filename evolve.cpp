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

//#define VERBOSE


// #define NOF_GENOMES 999
// #define GENE_VARIATION_PROB 0.5	//probability of variation when ancestor gene is aquired
// #define LOCUS_VARIATION_PROB 0.01 // 0.05  // probability of variating a nucleotide in a variated sequence
// #define GENE_DUPLICATION_PROB 0.001	//probability of duplicating a gene
// #define GENESET_VARIATION 0.01	//percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
// #define GENESET_VARIATION_ADD 0.9 //probability that the variation is a a gene add
// #define GENESET_VARIATION_REMOVE 0.1 //probability hat the variation if a gene removal
// #define RAND_SEED 123456789


void usage(std::string cmd){
    std::cout<<"Usage: "<<cmd<<" root_genome.peg hgt_pool.hgt oprefix tree.genome_parents sub_matrix GENE_VARIATION_PROB LOCUS_VARIATION_PROB GENE_DUPLICATION_PROB GENESET_VARIATION GENESET_VARIATION_ADD RAND_SEED\n";
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
                if(end >= g->sequence.size()){
                    end = g->sequence.size()-1;
                }
                //std::cout<<"--->"<<start<<" "<<end<<" "<<c<<"\n";
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
    };

    // std::vector<std::string>*
    // get_gene_sequences(){
    //     std::vector<std::string> *v = new std::vector<std::string>();
    //     for(auto & locus : this->loci){
    //         if(locus.strand == 1){
    //             std::string s = this->sequence.substr(locus.start, locus.end-locus.start);
    //             v->push_back(s);
    //         }
    //         else{
    //             std::string t = this->sequence.substr(locus.start, locus.end-locus.start);

    //             std::string s = this->sequence.substr(locus.start, locus.end-locus.start);
    //             for(int i=0; i<s.size(); i++){
    //                 s[i] = Genome::rc_symbol( t[ s.size()-1-i ] );
    //             }

    //             v->push_back(s);
    //         }
    //     }
    //     return v;
    // };


};


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

int int2codons_len = 58;
static const std::string int2codons[] = {
    "AAA",
    "AAC",
    "AAG",
    "AAT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "AGA",
    "AGC",
    "AGG",
    "AGT",
    "ATA",
    "ATC",
    //"ATG",
    "ATT",
    "CAA",
    "CAC",
    "CAG",
    "CAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "GTA",
    "GTC",
    //"GTG",
    "GTT",
    //"TAA",
    "TAC",
    //"TAG",
    "TAT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    //"TGA",
    "TGC",
    "TGG",
    "TGT",
    "TTA",
    "TTC",
    //"TTG",
    "TTT",
    };

bool to_revert(std::string &s, int p){
    //"ATG",
    //"GTG",
    //"TAA",
    //"TAG",
    //"TGA",
    //"TTG",
    if(s[p]=='A' && s[p+1]=='T' && s[p+2]=='G'){
        return true;
    }
    else if(s[p]=='G' && s[p+1]=='T' && s[p+2]=='G'){
        return true;
    }
    else if(s[p]=='T' && s[p+1]=='A' && s[p+2]=='A'){
        return true;
    }
    else if(s[p]=='T' && s[p+1]=='A' && s[p+2]=='G'){
        return true;
    }
    else if(s[p]=='T' && s[p+1]=='G' && s[p+2]=='A'){
        return true;
    }
    else if(s[p]=='T' && s[p+1]=='T' && s[p+2]=='G'){
        return true;
    }
    else{
        return false;
    }
};

static const char int2nuc[] = {'A','C','G','T'};

// std::string*
// make_random_sequence(int min_length, int max_length){

//     int length = randint(min_length, max_length) + 6;
//     char * c = new char[length + 1];
//     c[0] = 'A'; c[1] = 'T'; c[2] = 'G';
//     c[length-3] = 'T'; c[length-2] = 'A'; c[length-1] = 'A'; 
//     c[length] = 0;

//     for(int i=0; i<length; i++){
//         c[i] = int2nuc[ randint(0,4) ];
//     }

//     std::string *s = new std::string(c);
//     return s;
// };


std::string*
make_random_sequence(int min_length, int max_length){

    int length = randint(min_length, max_length) + 6;
    if(length%3==1) length+=2;
    if(length%3==2) length+=1;

    char * c = new char[length + 1];

    for(int i=0; i<length-3; i+=3){
        //c[i] = int2nuc[ randint(0,4) ];
        std::string s = int2codons[randint(0, int2codons_len)];
        
        c[i] = s[0];
        c[i+1] = s[1];
        c[i+2] = s[2];

    }

    c[0] = 'A'; c[1] = 'T'; c[2] = 'G';
    c[length-3] = 'T'; c[length-2] = 'A'; c[length-1] = 'A'; 
    c[length] = 0;

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


char** read_subsitution_matrix(std::string ifile){
    char ** m = new char*[4];
    for(int i=0; i<4; i++){
        m[i] = new char[100];
    }

    std::ifstream file(ifile);
    int x;

    for(int i=0; i<4; i++){
        int p=0;
        for(int k=0;k<100; k++){
            m[i][k] = int2nuc[i];
        }

        for(int j=0; j<4; j++){
            file >> x;
            //std::cout<<"@ "<<x<<" "<<int2nuc[j]<<"\n";
            for(int k=0;k<x && p<100; k++, p++){
                m[i][p] = int2nuc[j];
            }
        }
    }

    file.close();

    return m;
};
char substitute(char x, char** sub_matrix){
    if(x=='A') return sub_matrix[0][randint(0,99)];
    else if(x=='C') return sub_matrix[1][randint(0,99)];
    else if(x=='G') return sub_matrix[2][randint(0,99)];
    else if(x=='T') return sub_matrix[3][randint(0,99)];
    return 'A';
}


int main(int argc, char** argv){
    std::cout<<argc<<"\n";
    if(argc!=12){
        usage(argv[0]);
        return 0;
    }


    std::string itree = argv[4];
    std::string isunm = argv[5];
    float GENE_VARIATION_PROB = atof(argv[6]);
    float LOCUS_VARIATION_PROB = atof(argv[7]);
    float GENE_DUPLICATION_PROB = atof(argv[8]);
    float GENESET_VARIATION = atof(argv[9]);
    float GENESET_VARIATION_ADD = atof(argv[10]);
    float GENESET_VARIATION_REMOVE = 1.0 - GENESET_VARIATION_ADD;
    int RAND_SEED = atoi(argv[11]);

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
        //std::cout<<"@@@ "<<l.start<<" "<<l.end<<" "<<length<<" "<<(length %3)<<"\n";
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
    std::cout<<"Reading population tree...\n";
    std::cout<<argv[4]<<"\n";
    std::map<int, std::vector<int> > tchilds;
    std::map<int,int> genome_parents;
    std::ifstream tfile(argv[4]);
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
    std::cout<<argv[5]<<"\n";
    char **subm = read_subsitution_matrix(argv[5]);

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

    //std::vector<Genome*> genomes;
    std::map<int, Genome*> genomes;

    //genomes.push_back(root_genome);
    genomes[0] = root_genome;
    //genome_parents[0] = -1;
    for(Locus &l : root_genome->loci){
        gene_parents[ std::pair<int,int>(0, l.id) ] = std::pair<int,int>(-1, -1);
    }


    std::mt19937_64 rng;
    rng.seed(RAND_SEED);
    std::uniform_real_distribution<double> unif(0, 1);

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
            if( unif(rng) <= GENE_VARIATION_PROB){
                
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
                        if( unif(rng) <= LOCUS_VARIATION_PROB){
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
            if( unif(rng) <= GENE_DUPLICATION_PROB ){
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
        std::cout<<"gene variation set size is "<<std::ceil(  parent_genome->loci.size() * GENESET_VARIATION )<<" / "<<parent_genome->loci.size()<<"\n";
#endif
        for(int gv = 0; gv<std::ceil(  parent_genome->loci.size() * GENESET_VARIATION ); gv++){
            if( unif(rng) <= GENESET_VARIATION_REMOVE){
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
                                new_genome->sequence.substr(0, range_start) + 
                                new_genome->sequence.substr(i);

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
                        new_genome->sequence.substr(0, range_start) + 
                        new_genome->sequence.substr(i);

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
            if( unif(rng) <= GENESET_VARIATION_ADD){
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

        //genomes.push_back(new_genome);
        genomes[current_genome_id] = new_genome;

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

};
