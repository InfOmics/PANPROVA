#ifndef GENE_UTIL_HH
#define GENE_UTIL_HH

#include <string>
#include <vector>
#include "codon.hh"
#include "randoms.hh"

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

#endif