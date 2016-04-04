//
// Created by klimas on 04.04.16.
//

#ifndef BIO2_RANDOMHEURISTIC_H
#define BIO2_RANDOMHEURISTIC_H


#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <random>
#include <array>
#include <algorithm>



#define A_BIN 0
#define C_BIN 1
#define G_BIN 2
#define T_BIN 3

#define OLI_LENGTH 10
#define MATRIX_COUNT 1048576
#define MAX_NEGATIVE 9

//#define PRINT_RESULT
#define PRINT_STATS

using namespace std;

class RandomHeuristic{
private:
    unsigned m; // number of nucleotides in file
    array<int,MAX_NEGATIVE> countOffset = {{0}};
    unsigned negativeError = 0;
    deque<unsigned> result;
    bool microArray[MATRIX_COUNT] = {false};

    //randomVector
    array< vector<int>,MAX_NEGATIVE> randomVectors;
    vector<unsigned > randomVector;
    unsigned randomVectorIterator;
    array<unsigned,MAX_NEGATIVE> randomVectorLimit;

public:
    void run();

    void runRandom();
    int randomSuccessor(uint32_t base, unsigned offset);
    int randomPredecessor(uint32_t nucleotide, unsigned offset);

    void initRandomVector();
    void shuffleVector(int offset);
    void generateRandomVector(unsigned min, unsigned max);


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


#endif //BIO2_RANDOMHEURISTIC_H
