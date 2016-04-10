//
// Created by klimas on 05.04.16.
//

#include "MonteCarlo.h"


int MonteCarlo::successor(uint32_t nucleotide, unsigned offset) {
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

int MonteCarlo::predecessor(uint32_t nucleotide, unsigned offset){
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


int MonteCarlo::recursiveSuccessor(uint32_t nucleotide, unsigned maxOffset) {
    int cRate = 0;
    for (int i = 1; i <= maxOffset; ++i) {
        int base = (nucleotide << (2 * i)) & ((1 << 20) - 1);
        for (int j = 0; j < pow(4, i); ++j) {
            int r = base + randomVectors[i - 1].at(j);
            if (microArray[r]) {
                result.push_back(r);
                microArray[r] = false;
                negativeError+=i-1;
                if ((++countOffset[i - 1]) == randomVectorLimit[i - 1]) {
                    shuffleVector(i);
                    randomVectorIterators[i-1] = randomVectorIterators[i-1] >= randomVector.size() ? 0 : randomVectorIterators[i-1];
                    randomVectorLimit[i - 1] += randomVector[randomVectorIterators[i-1]++];
                }
                cRate = recursiveSuccessor(r, maxOffset);
                --countOffset[i - 1];
                negativeError-=i-1;
                microArray[r] = true;
                result.pop_back();
                return cRate;
            }
        }
    }
//    printResultAsString();
//    printStats();
    return rate();
}

void MonteCarlo::runForward() {
    cout << "wynik: " << endl;
//    cout << recursiveSuccessor(result.at(0), 2) << endl; //musi być mniejszy niż MAX_NEGATIVE
    cout << runN(result.at(0), 2,50,[this](uint32_t nucleotide, unsigned offset){ return recursiveSuccessor(nucleotide,offset);}) << endl; //musi być mniejszy niż MAX_NEGATIVE

}

