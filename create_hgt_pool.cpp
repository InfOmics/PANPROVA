#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cctype>
#include <stdlib.h>
#include <iterator>



#define JACCARD_K 6
#define JACCARD_ROOT_SIM 0.3
#define JACCARD_HGT_SIM 0.5


void usage(std::string cmd){
    std::cout<<"Usage: "<<cmd<<" ilistfile ofile\n";
    std::cout<<"ilistifle contains paths to PEG files, one per line. The first one is the root genome.";
}



class Locus{
public:
    int start;
    int end;
    int strand;// 1, -1

    Locus(int _start, int _end, int _strand)
        : start(_start), end(_end), strand(_strand)
    {}

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
        return os << "("<<tc.start<<","<<tc.end<<","<<tc.strand<<")";
    }
};


class Genome{
public:
    std::string name;
    std::string sequence;
    std::vector<Locus> loci;

    Genome(){
        name = "";
        sequence = "";
    }

    static 
    Genome* read_from_file(std::string &ifile) {
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

        int c;
        int state = 0;
        int start, end;
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
                g->loci.push_back( Locus(start,end,c) );
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
                std::string s = this->sequence.substr(locus.start, locus.end-locus.start);
                v->push_back(s);
            }
            else{
                std::string t = this->sequence.substr(locus.start, locus.end-locus.start);

                std::string s = this->sequence.substr(locus.start, locus.end-locus.start);
                for(int i=0; i<s.size(); i++){
                    s[i] = Genome::rc_symbol( t[ s.size()-1-i ] );
                }

                v->push_back(s);
            }
        }
        return v;
    };


};


std::map<std::string, int> *
get_kmer_multiplicities(std::string &s, int k){
    int *N = new int[s.size()];
    for(int i=0; i<s.size(); i++){
        N[i] = 0;
    }
    int lastN = s.size();
    for(int i=s.size()-1; i>=0; i--){
        if((s[i]!='A') & (s[i]!='T') & (s[i]!='C') & (s[i]!='G')){
            lastN = i;
        }
        N[i] = lastN - i;
    }

    std::map<std::string, int> * mults = new std::map<std::string, int>();

    for(int i=0; i<=s.size()-k+1; i++){
        if(N[i] >= k){
            std::string kmer = s.substr(i,k);
            (*mults)[kmer] = (*mults)[kmer] + 1;
        }
    }

    delete [] N;

    return mults;
};

double
generalized_jaccard(std::map<std::string, int> &v1, std::map<std::string, int> &v2){
    double num = 0.0;
    double den = 0.0;

    int v2visit = 0;
    int c2;
    
    for(const auto &p : v1){
        c2 = v2[p.first];
        v2visit += c2;

        if(p.second >= c2){
            den += p.second;
            num += c2;
        }
        else{
            num += p.second;
            den += c2;
        }
    }
    
    int v2sum = 0;
    for(const auto &p : v2){
        v2sum += p.second;
    }


    return num / (den + (v2sum - v2visit));
}



int main(int argc, char** argv){
    std::cout<<argc<<"\n";
    if(argc!=3){
        usage(argv[0]);
        return 0;
    }

    //int jaccard_k = 6;
    //double jaccard_root_sim = 0.3;
    //double jaccard_hgt_sim = 0.5;


    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading root genome...\n";
    std::ifstream file(argv[1]);
    std::string line; 
    std::getline(file, line);
    std::cout<<line<<"\n";
    Genome *root_genome = Genome::read_from_file(line);
    file.close();

    std::cout<<"Root genome is "<<root_genome->sequence.size()<<" nucleotides long, with "<<root_genome->loci.size()<<" genetic loci\n";

    std::cout<<"----------------------------------------\n";
    std::cout<<"Reading HGT genomes...\n";
    std::vector<Genome*> hgt_genomes;
    file.open(argv[1]);
    std::getline(file, line);
    while(std::getline(file, line)){
        std::cout<<line<<"\n";
        Genome *hgt = Genome::read_from_file(line);
        if(hgt){
            std::cout<<hgt->sequence.size()<<" "<<hgt->loci.size()
                     <<" "<<*(hgt->loci.begin())<<" \n";
            std::cout<<hgt->sequence.substr(0,100)<<"\n";
            hgt->name = line;
            hgt_genomes.push_back(hgt);
        }else{
            std::cout<<"genome not loaded.\n";
        }
        std::cout<<"--------------------\n";
    }
    file.close();
    std::cout<<"----------------------------------------\n";

    std::cout<<"Retrieving gene sequences ...\n";
    std::cout<<"--------------------\n";
    std::vector<std::string> *root_genes = root_genome->get_gene_sequences();
    std::cout<<"retrieved "<<root_genes->size()<<" root genes\n";
    std::cout<<"--------------------\n";

    std::vector<std::string> genome_names;
    std::vector<int> genes_to_genome;
    std::vector<Locus> loci;

    std::vector<std::string> *hgt_genes = new std::vector<std::string>();
    int gid = 0;
    for(auto &genome : hgt_genomes){
        genome_names.push_back(genome->name);
        std::vector<std::string> *genes = genome->get_gene_sequences();
        hgt_genes->insert(hgt_genes->end(), genes->begin(), genes->end());
        for(int i=0; i<hgt_genes->size(); i++){
            genes_to_genome.push_back(gid);
        }
        for(auto & locus : genome->loci){
            loci.push_back(locus);
        }
        gid++;
    }
    std::cout<<"retrieved an initial pool of "<<hgt_genes->size()<<" HGT genes\n";

    delete root_genome;
    for(auto & e : hgt_genomes){
        delete e;
    }

    std::cout<<"----------------------------------------\n";

    std::cout<<"Discarding HGT genes by their similarity with root genes ...\n";
    std::cout<<"--------------------\n";


    std::cout<<"indexing root genes ...\n";
    
    std::map<std::string, double*> root_kmer_mults;
    double *root_kmer_sums = new double[root_genes->size()];

    int genei = 0;
    for(auto &s : *root_genes){
        std::map<std::string, int> * mults = get_kmer_multiplicities(s,JACCARD_K);
        root_kmer_sums[genei] = 0;
        for(const auto &p : *mults){
//            std::cout<<p.first<<" "<<p.second<<"\n";
            if(root_kmer_mults.find(p.first) == root_kmer_mults.end()){
                //int *x = new int[ root_genes->size()];
                //x[genei] = p.second;

                root_kmer_mults[ p.first ] =  (double*)calloc( root_genes->size(), sizeof(double));

                /*
                root_kmer_mults[ p.first ] = new int[ root_genes->size()];
                for(int i=0; i<root_genes->size(); i++){
                    root_kmer_mults[ p.first ][i] = 0;
                }
                */
            }
            //else{
                root_kmer_mults[ p.first ][genei] = p.second;
            //}
            root_kmer_sums[genei] += p.second;

            //std::cout<<genei<<" "<<p.first<<" "<<p.second<<" "<<root_kmer_mults[ p.first ][genei] <<"\n";
        }
        genei++;
    }

    
    std::cout<<root_kmer_mults.size()<<" "<<JACCARD_K<<"-mers found\n";
    /*std::cout<<"--------------------\n";

    for(const auto &p : root_kmer_mults){
        std::cout<<p.first<<"\n";
        for(int i=0; i<10; i++){
            std::cout<<p.second[i]<<" ";
        }
        std::cout<<"\n";
    }*/

    
    std::cout<<"--------------------\n";

    double * nums = new double[root_genes->size()];
    double * dens = new double[root_genes->size()];
    double * sums = new double[root_genes->size()];
    double *root_kmers;
    std::map<std::string, double*>::iterator root_kmers_it;

    std::map<std::string, double*> hgt_kmer_mults;
    double *hgt_kmer_sums = new double[hgt_genes->size()];
    for(int i=0; i<hgt_genes->size(); i++){
        hgt_kmer_sums[i] = 0;
    }

    std::set<int> hgt_selection;
    bool deleteit;
    double jvalue;

    std::cout<<"comparing hgt genes to root genes ...\n";
    for(int h=0; h<hgt_genes->size(); h++){
        std::cout<<h<<"/"<<hgt_genes->size()<<"\n";
        std::map<std::string, int> * mults = get_kmer_multiplicities((*hgt_genes)[h],JACCARD_K);
        
        for(int i=0; i<root_genes->size(); i++){
            nums[i] = dens[i] = sums[i] = 0;
        }

    
        for(const auto &p : *mults){
            //if(root_kmer_mults.find(p.first) != root_kmer_mults.end()){
                //root_kmers = root_kmer_mults[p.first];
            root_kmers_it = root_kmer_mults.find(p.first);
            if(root_kmers_it != root_kmer_mults.end()){
                //std::cout<<p.first<<"\n";
                root_kmers = root_kmers_it->second;
                
                //for(int i=0; i<10; i++){
                //    std::cout<<root_kmers[i]<<" ";
                //}
                //std::cout<<"\n";

                for(int i=0; i<root_genes->size(); i++){
                    if(p.second > root_kmers[i]){
                        nums[i] += root_kmers[i];
                        dens[i] += p.second;
                    }
                    else{
                        dens[i] += root_kmers[i];
                        nums[i] += p.second;
                    }
                    sums[i] += root_kmers[i];
                }
            }
            else{
                for(int i=0; i<root_genes->size(); i++){
                    dens[i] += p.second;
                }
            }
        }

        deleteit = false;
        for(int i=0; i<root_genes->size(); i++){
            jvalue = nums[i] / ( dens[i] + (root_kmer_sums[i] - sums[i]) );
            if(jvalue >= JACCARD_ROOT_SIM){
                deleteit = true;
                break;
            }
        }

        if(!deleteit){
            hgt_selection.insert(h);

            for(const auto &p : *mults){
                if(hgt_kmer_mults.find(p.first) == hgt_kmer_mults.end()){
                    hgt_kmer_mults[ p.first ] = (double*)calloc( hgt_genes->size(), sizeof(double));
                    /*
                    hgt_kmer_mults[ p.first ] = new int[ hgt_genes->size()];
                    for(int i=0; i<hgt_genes->size(); i++){
                        hgt_kmer_mults[ p.first ][i] = 0;
                    }
                    */
                }
                hgt_kmer_mults[ p.first ][h] = p.second;
                hgt_kmer_sums[h] += p.second;
            }
        }
    }


    delete [] nums;
    delete [] dens;
    delete [] sums;
    for(auto &p : root_kmer_mults){
        free(p.second);
    }

    std::cout<<"selected "<<hgt_selection.size()<<"/"<< hgt_genes->size() <<" valid HGT genes\n";
    std::cout<<"----------------------------------------\n";
    std::cout<<"Discarding HGT genes by similarity between them ...\n";
    std::cout<<"--------------------\n";



    std::cout<<"calculating jaccard similarities ...\n";
    double **nums_map = new double*[hgt_genes->size()];
    double **dens_map = new double*[hgt_genes->size()];
    for(int i=0; i<hgt_genes->size(); i++){
        nums_map[i] = (double*)calloc(hgt_genes->size(), sizeof(double)); //new double[hgt_genes->size()];
        dens_map[i] = (double*)calloc(hgt_genes->size(), sizeof(double)); //new double[hgt_genes->size()];
    }


    int hgt_sel_length = hgt_selection.size();
    int * hgt_sel = new int[hgt_selection.size()];
    int h1 = 0, h2, h1i, h2i;
    for(auto &h : hgt_selection){
        hgt_sel[h1] = h;
        h1++;
    }
    std::pair<int,int> hh(0,0);
    double hs, hsmax;

    for(const auto &p : hgt_kmer_mults){
        std::cout<<p.first<<"\n";
        for(h1i=0; h1i<hgt_sel_length; h1i++){
            h1 = hgt_sel[h1i];
            if( (p.second[h1]>0)) {
                //std::cout<<p.first<<" "<<h1i<<" "<<h1<<"\n";

                for(h2i=h1i+1; h2i<hgt_sel_length; h2i++){
                    h2 = hgt_sel[h2i];

                    if(p.second[h1] > p.second[h2]){
                        nums_map[h1][h2] += p.second[h2];
                        dens_map[h1][h2] += p.second[h1];
                    }
                    else{
                        nums_map[h1][h2] += p.second[h1];
                        dens_map[h1][h2] += p.second[h2];
                    }
                }
            }
        }
    }


    std::cout<<"discarding similar genes ...\n";
    std::set<int> todelete;

    for(h1i=0; h1i<hgt_sel_length; h1i++){
        h1 = hgt_sel[h1i];
        for(h2i=h1i+1; h2i<hgt_sel_length; h2i++){
            h2 = hgt_sel[h2i];

            if(dens_map[h1][h2] > 0.0){
                jvalue = nums_map[h1][h2] / dens_map[h1][h2];
                if(jvalue >= JACCARD_HGT_SIM){
                    todelete.insert(h2);
                }
            }
        }
    }

    std::cout<<todelete.size()<< " HGT genes discarded \n";

    std::cout<<"----------------------------------------\n";

    std::cout<<"Writing HGT genes to "<<argv[2]<<" and "<<argv[2]<<".fa\n";

    std::ofstream myfile;
    myfile.open (argv[2]);

    std::ofstream myfile2;
    myfile2.open (std::string(argv[2])+".fa");


    //std::vector<std::string> towrite;
    int count = 0;
    for(h1i=0; h1i<hgt_sel_length; h1i++){
        h1 = hgt_sel[h1i];
        if(todelete.find(h1) == todelete.end()){
            //towrite.push_back(hgt_genes[h1]);
            myfile << (*hgt_genes)[h1]<<"\n";

            myfile2 <<">"<< genome_names[genes_to_genome[h1]]<<":"<< loci[h1].start<<":"<<loci[1].end<<":"<<loci[h1].strand  <<"\n";
            myfile2 << (*hgt_genes)[h1]<<"\n";

            count++;
        }
    }
    myfile.flush();
    myfile.close();
    myfile2.flush();
    myfile2.close();
    std::cout<<count<<" genes written on file\n";
    std::cout<<"----------------------------------------\n";
    std::cout<<"done\n";


    std::cout<<"----------------------------------------\n";

    return 0;
}