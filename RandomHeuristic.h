//
// Created by klimas on 04.04.16.
//

#ifndef BIO2_RANDOMHEURISTIC_H
#define BIO2_RANDOMHEURISTIC_H


#include "Heuristic.h"

class RandomHeuristic:public Heuristic{

public:

    int successor(uint32_t base, unsigned offset) override final ;
    int predecessor(uint32_t nucleotide, unsigned offset) override final ;

};


#endif //BIO2_RANDOMHEURISTIC_H
