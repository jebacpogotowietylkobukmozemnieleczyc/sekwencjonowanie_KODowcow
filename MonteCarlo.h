//
// Created by klimas on 05.04.16.
//

#ifndef BIO2_MONTECARLO_H
#define BIO2_MONTECARLO_H


#include "Heuristic.h"

class MonteCarlo : public Heuristic {
public:
    //todo monte carlo
    void run() override final ;
    int successor(uint32_t nucleotide, unsigned offset);
    int predecessor(uint32_t nucleotide, unsigned offset);

};


#endif //BIO2_MONTECARLO_H
