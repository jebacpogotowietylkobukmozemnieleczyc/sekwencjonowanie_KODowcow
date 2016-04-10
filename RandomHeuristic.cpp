//
// Created by klimas on 04.04.16.
//

#include "RandomHeuristic.h"


int RandomHeuristic::successor(uint32_t nucleotide, unsigned offset) {
    int base = (nucleotide << (2 * offset)) & ((1 << 20) - 1);
    for (int j = 0; j <  pow(4,offset); ++j) {
        int r = base + randomVectors[offset - 1].at(j);
        if(microArray[r]) {
            result.push_back(r);
            microArray[r] = false;
            return r;
        }
    }
    return -1;
}

int RandomHeuristic::predecessor(uint32_t nucleotide, unsigned offset){
    int step = 262144/pow(4,offset-1);
    int base = (nucleotide>>(2*offset)) & ((1<<20)-1);
    for (int j = 0; j <  pow(4,offset); ++j) {
        int r = base + randomVectors[offset - 1].at(j) * step;
        if(microArray[r]) {
            microArray[r] = false;
            result.push_front(r);
            return r;
        }
    }

    return -1;
}




