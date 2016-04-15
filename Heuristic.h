//
// Created by klimas on 04.04.16.
//

#ifndef BIO2_HEURISTIC_H
#define BIO2_HEURISTIC_H


#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <random>
#include <array>
#include <algorithm>
#include <functional>
#include <memory>


#define A_BIN 0
#define C_BIN 1
#define G_BIN 2
#define T_BIN 3

#define OLI_LENGTH 10
#define MATRIX_COUNT 1048576
#define MAX_NEGATIVE 9

#define PRINT_RESULT
#define PRINT_STATS

using namespace std;


class Heuristic{
protected:
    unsigned m = 0 ; // number of nucleotides in file
    //todo zero?
    array<int,MAX_NEGATIVE> countOffset;
    unsigned negativeError = 0;
    deque<unsigned> result;

protected:
    array<bool,MATRIX_COUNT> microArray;

    //randomVector
    array< vector<int>,MAX_NEGATIVE> randomVectors;
    vector<unsigned > randomVector;
    array<unsigned,MAX_NEGATIVE> randomVectorLimit;
    array<unsigned,MAX_NEGATIVE> randomVectorIterators;

public:
    Heuristic() {
        countOffset.fill(0);
        microArray.fill(false);
        randomVectorLimit.fill(0);
    }

    int runN(uint32_t nucleotide, unsigned maxOffset, unsigned n , function <int(uint32_t, unsigned)> f);
    void run();

    virtual int successor(uint32_t base, unsigned offset) = 0;
    virtual int predecessor(uint32_t nucleotide, unsigned offset) = 0;

    void initRandomVector();
    void shuffleVector(int offset);
    void generateRandomVector(unsigned min, unsigned max);
    void initRandomVectorLimit();


    unsigned rate();
    bool checkIfMatch(uint32_t leftNucleotide, uint32_t rightNucleotide);
    bool checkIfMatch(uint32_t leftNucleotide, uint32_t rightNucleotide, uint32_t offset);
    uint32_t getRandom(uint32_t min, uint32_t max);
    bool test();
    void printStats();
    uint32_t stringIntoIntCoder(string sequence);
    string intIntoStringCoder(uint32_t codedNumber);
    void readFromFile(const char *fname, unsigned fileLength);
    void printResult();
    void printResultAsString();


};


#endif //BIO2_HEURISTIC_H
