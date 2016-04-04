//
// Created by klimas on 04.04.16.
//

#include "RandomHeuristic.h"


void RandomHeuristic::run() {
    if(result.size()!=1){
        cout << "cos sie zepsulo" << std::endl;
    }
    uint32_t successorNucleotide = result.at(0);
    uint32_t predecessorNucleotide = result.at(0);



    int nextNucleotide = -1;
    for (int i = 1; i <= MAX_NEGATIVE; ++i) {
        while ((nextNucleotide = randomSuccessor(successorNucleotide, i)) != -1) {
            if ((++countOffset[i - 1]) == randomVectorLimit[i - 1]) {
                shuffleVector(i);
                randomVectorIterator = randomVectorIterator >= randomVector.size() ? 0 : randomVectorIterator;
                randomVectorLimit[i - 1] += randomVector[randomVectorIterator++];
            }
            if (i > 1){
                --i;
                negativeError+=i-1;
            }
            successorNucleotide = nextNucleotide;
        };


        while ((nextNucleotide = randomPredecessor(predecessorNucleotide, i)) != -1) {
            if ((++countOffset[i - 1]) == randomVectorLimit[i - 1]) {
                shuffleVector(i);
                randomVectorIterator = randomVectorIterator >= randomVector.size() ? 0 : randomVectorIterator;
                randomVectorLimit[i - 1] += randomVector[randomVectorIterator++];
            }
            if (i > 1){
                --i;
                negativeError+=i-1;
            }
            predecessorNucleotide = nextNucleotide;
        }

    }

}

int RandomHeuristic::randomSuccessor(uint32_t nucleotide, unsigned offset) {
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

int RandomHeuristic::randomPredecessor(uint32_t nucleotide, unsigned offset){
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


void RandomHeuristic::runForward() {
    cout << "wynik: " << endl;
    cout << recursiveSuccessor(result.at(0), 2) << endl; //musi być mniejszy niż MAX_NEGATIVE
//    cout << runN(result.at(0), 2,50,recursiveSuccessor) << endl; //musi być mniejszy niż MAX_NEGATIVE

}


int RandomHeuristic::recursiveSuccessor(uint32_t nucleotide, unsigned maxOffset) {
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
                    randomVectorIterator = randomVectorIterator >= randomVector.size() ? 0 : randomVectorIterator;
                    randomVectorLimit[i - 1] += randomVector[randomVectorIterator++];

                    cRate = recursiveSuccessor(r, maxOffset);

                    randomVectorIterator = randomVectorIterator <= 0 ? 0 : randomVectorIterator-1;
                    randomVectorLimit[i - 1] -= randomVector[randomVectorIterator++];
                    ++countOffset[i - 1];
                }
                else{
                    cRate = recursiveSuccessor(r, maxOffset);
                }
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


int RandomHeuristic::runN(uint32_t nucleotide, unsigned maxOffset, unsigned n, function<int(uint32_t, unsigned)> f) {
    int maxRate;
    for (int i = 0; i < n; ++i) {
        maxRate = max(maxRate, f(nucleotide, maxOffset));
    }
    return maxRate;
}
