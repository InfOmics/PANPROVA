#ifndef GENOME_HH
#define GENOME_HH

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "Locus.hh"

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


#endif