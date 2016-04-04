//
// Created by klimas on 04.04.16.
//

#ifndef BIO2_RANDOMHEURISTIC_H
#define BIO2_RANDOMHEURISTIC_H


#include "Heuristic.h"

class RandomHeuristic:public Heuristic{

public:
    void run() override final;

    void runForward();
    int randomSuccessor(uint32_t base, unsigned offset);
    int randomPredecessor(uint32_t nucleotide, unsigned offset);


    int runN(uint32_t nucleotide, unsigned maxOffset, unsigned n ,function <int(uint32_t,unsigned)> f);
    int recursiveSuccessor(uint32_t nucleotide, unsigned maxOffset);
    //todo recursicePredecessor
    int recursicePredecessor(uint32_t nucleotide, unsigned maxOffset);

};


#endif //BIO2_RANDOMHEURISTIC_H
