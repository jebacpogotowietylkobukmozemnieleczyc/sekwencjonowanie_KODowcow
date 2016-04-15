// // Created by klimas on 05.04.16.
//

#include "MonteCarlo.h"


int MonteCarlo::successor(uint32_t nucleotide, unsigned offset) {
    tuple<uint32_t, unsigned> maxRate{0, 0};
    vector<uint32_t> successors;
    int base = (nucleotide << (2 * offset)) & ((1 << 20) - 1);
    for (int j = base; j < pow(4, offset); ++j) {
        if (microArray[j]) {
            successors.push_back(j);
        }
    }
    if (successors.size() > 1) {
        for (uint32_t element :successors) {
            tuple<uint32_t, unsigned> newRate{element, runForward(element)};
            maxRate = max(maxRate, newRate,
                          [](const tuple<uint32_t, unsigned> &tuple1, const tuple<uint32_t, unsigned> &tuple2) {
                              return get<1>(tuple1) < get<1>(tuple2);
                          });
        }
        result.push_back(get<0>(maxRate));
        microArray[get<0>(maxRate)] = false;
        return get<0>(maxRate);
    }
    else if (successors.size() == 1) {
        uint32_t newNucleotide = successors.at(0);
        result.push_back(newNucleotide);
        microArray[newNucleotide] = false;
        return newNucleotide;
    }
    return -1;
}

int MonteCarlo::predecessor(uint32_t nucleotide, unsigned offset) {
    tuple<uint32_t, unsigned> maxRate{0, 0};
    vector<uint32_t> predecessors;

    int step = 262144 / pow(4, offset - 1);
    int base = (nucleotide >> (2 * offset)) & ((1 << 20) - 1);
    for (int j = base; j < MATRIX_COUNT; j += step) {
        if (microArray[j]) {
            predecessors.push_back(j);
        }
    }

    if (predecessors.size() > 1) {
        for (uint32_t element :predecessors) {
            tuple<uint32_t, unsigned> newRate{element, runBackward(element)};
            maxRate = max(maxRate, newRate,
                          [](const tuple<uint32_t, unsigned> &tuple1, const tuple<uint32_t, unsigned> &tuple2) {
                              return get<1>(tuple1) < get<1>(tuple2);
                          });
        }
        result.push_front(get<0>(maxRate));
        microArray[get<0>(maxRate)] = false;
        return get<0>(maxRate);
    }
    else if (predecessors.size() == 1) {
        uint32_t newNucleotide = predecessors.at(0);
        result.push_front(newNucleotide);
        microArray[newNucleotide] = false;
        return newNucleotide;
    }
    return -1;
}


int MonteCarlo::recursivePredecessor(uint32_t nucleotide, unsigned maxOffset){
    int cRate = 0;
    int step = 262144;
    for (int i = 1; i <= maxOffset; ++i) {
        int base = (nucleotide >> (2 * i)) & ((1 << 20) - 1);
        for (int j = 0; j < pow(4, i); ++j) {
            int r = base + randomVectors[i - 1].at(j) * step;
            if (microArray[r]) {
                result.push_front(r);
                microArray[r] = false;
                negativeError += i - 1;
                if ((++countOffset[i - 1]) == randomVectorLimit[i - 1]) {
                    shuffleVector(i);
                    randomVectorIterators[i - 1] =
                            randomVectorIterators[i - 1] >= randomVector.size() ? 0 : randomVectorIterators[i - 1];
                    randomVectorLimit[i - 1] += randomVector[randomVectorIterators[i - 1]++];
                }
                cRate = recursiveSuccessor(r, maxOffset);
                --countOffset[i - 1];
                negativeError -= i - 1;
                microArray[r] = true;
                result.pop_front();
                return cRate;
            }
        }
        step /= 4;
    }
//    printResultAsString();
//    printStats();
    return rate();
}

int MonteCarlo::recursiveSuccessor(uint32_t nucleotide, unsigned maxOffset){
    int cRate = 0;
    for (int i = 1; i <= maxOffset; ++i) {
        int base = (nucleotide << (2 * i)) & ((1 << 20) - 1);
        for (int j = 0; j < pow(4, i); ++j) {
            int r = base + randomVectors[i - 1].at(j);
            if (microArray[r]) {
                result.push_back(r);
                microArray[r] = false;
                negativeError += i - 1;
                if ((++countOffset[i - 1]) == randomVectorLimit[i - 1]) {
                    shuffleVector(i);
                    randomVectorIterators[i - 1] =
                            randomVectorIterators[i - 1] >= randomVector.size() ? 0 : randomVectorIterators[i - 1];
                    randomVectorLimit[i - 1] += randomVector[randomVectorIterators[i - 1]++];
                }
                cRate = recursiveSuccessor(r, maxOffset);
                --countOffset[i - 1];
                negativeError -= i - 1;
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

unsigned MonteCarlo::runForward(uint32_t nucleotide) {
    //todo //musi być mniejszy niż MAX_NEGATIVE
    return runN(nucleotide, 2, 50, [this](uint32_t nucleotide, unsigned offset) {
        return recursiveSuccessor(nucleotide, offset);
    }); //musi być mniejszy niż MAX_NEGATIVE

}

unsigned MonteCarlo::runBackward(uint32_t nucleotide) {
    //todo //musi być mniejszy niż MAX_NEGATIVE
    return runN(nucleotide, 2, 50, [this](uint32_t nucleotide, unsigned offset) {
        return recursivePredecessor(nucleotide, offset);
    }); //musi być mniejszy niż MAX_NEGATIVE

}
