//
// Created by klimas on 05.04.16.
//

#ifndef BIO2_MONTECARLO_H
#define BIO2_MONTECARLO_H


#include "Heuristic.h"

class MonteCarlo : public Heuristic {
public:
    int successor(uint32_t nucleotide, unsigned offset) override final;
    int predecessor(uint32_t nucleotide, unsigned offset) override final;


    int recursiveSuccessor(uint32_t nucleotide, unsigned maxOffset);
    //todo recursicePredecessor
    int recursicePredecessor(uint32_t nucleotide, unsigned maxOffset);

    void runForward();

};


#endif //BIO2_MONTECARLO_H
