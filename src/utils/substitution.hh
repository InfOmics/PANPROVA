#ifndef SUBSTITUTION_UTIL_HH
#define SUBSTITUTION_UTIL_HH

#include <string>
#include <fstream>
#include "randoms.hh"

static const char int2nuc[] = {'A','C','G','T'};

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

#endif