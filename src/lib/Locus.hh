#ifndef LOCUS_HH
#define LOCUS_HH

#include <iostream>
#include <fstream>

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


#endif