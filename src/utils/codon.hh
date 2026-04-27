#ifndef CODON_UTIL_HH
#define CODON_UTIL_HH

#include <string>


// 64 all possible codons except the 6 special ones, so 64-6=58:
// the three start codons (ATG, GTG, TTG) and the three stop codons
// (TAA, TAG, TGA)
const int int2codons_len = 58;
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


#endif