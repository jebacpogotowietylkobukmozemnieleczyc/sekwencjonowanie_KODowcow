#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <random>

#include "timer.hpp"

#define A_BIN 0
#define C_BIN 1
#define G_BIN 2
#define T_BIN 3

#define OLI_LENGTH 10
#define MATRIX_COUNT 1048576
#define MAX_NEGATIVE 1

//#define PRINT_RESULT
#define PRINT_STATS
#define COUNT_STATS

using namespace std;

uint32_t getRandom(uint32_t min, uint32_t max){
    std::uniform_int_distribution<mt19937::result_type > udist(min, max);
    mt19937 rand;
    random_device randomDevice;
    mt19937::result_type const seedval = randomDevice(); // get this from somewhere
    rand.seed(seedval);

    mt19937::result_type random_number = udist(rand);

    return random_number;

}


class OlinukleoLibrary{
private:
    vector<int> oliLib;
    unsigned m;

public:
    unsigned negativeOffset[MAX_NEGATIVE] = {0};
    unsigned negativeError = 0;
    deque<int> result;
    bool microArray[MATRIX_COUNT] = {false};
    uint32_t stringIntoIntCoder(string sequence);
    string intIntoStringCoder(uint32_t codedNumber);
    int readFromFile(const char* fname, unsigned fileLength);
    void printVector();
    void printVectorAsString();

    void greedyAlgorithm(int startNukleo);
    int successor(uint32_t base, unsigned offset);
    int predecessor(uint32_t nucleotide, unsigned offset);
    bool test();


    void printStats();
};

uint32_t OlinukleoLibrary::stringIntoIntCoder(string sequence){
    uint32_t result = 0;

    for( int i = 0; i < sequence.length(); i++){
        switch(sequence.at(i)){
            case 'A': result = result | A_BIN;
                break;

            case 'C': result = result | C_BIN;
                break;

            case 'G': result = result | G_BIN;
                break;

            case 'T': result = result | T_BIN;
                break;
        }
        result = result<<2;
    }
    result = result>>2;

    return result;
}


string OlinukleoLibrary::intIntoStringCoder(uint32_t codedNumber){
    string result;

    for( int i = OLI_LENGTH-1; i >= 0; i--){
        switch( (codedNumber>>(i*2))&3){
            case A_BIN: result+="A";
                break;

            case C_BIN: result+="C";
                break;

            case G_BIN: result+="G";
                break;

            case T_BIN: result+="T";
                break;
        }
    }

    return result;
}

int OlinukleoLibrary::readFromFile(const char* fname, unsigned fileLength){
    fstream file;
    int iterator = 0;
    string tmp;
    int randomNukleo = -1;

    int randomLine = getRandom(0, fileLength - 1);
    file.open(fname, ios::in);
    while (!file.eof()){
        file >> tmp;
        microArray[OlinukleoLibrary::stringIntoIntCoder(tmp)]=true;
        //todo take properties from file name
        if(iterator==randomLine){
            randomNukleo=OlinukleoLibrary::stringIntoIntCoder(tmp);
            microArray[OlinukleoLibrary::stringIntoIntCoder(tmp)]= false;
        }
        iterator++;
    }

    file.close();

    //todo empty line at end
    m = iterator - 1;
    return randomNukleo;
}

void OlinukleoLibrary::printVector(){
    for (int i = 0; i < oliLib.size(); i++){
        cout << oliLib[i] << endl;
    }
}

void OlinukleoLibrary::printVectorAsString(){
    for (int i = 0; i < oliLib.size(); i++){
        cout << OlinukleoLibrary::intIntoStringCoder(oliLib[i]) << endl;
    }
}


bool checkIfMatch(uint32_t leftOli, uint32_t rightOli){
    return (leftOli & ((1 << 18) - 1) ) == ((rightOli >> 2) & ((1 << 18) - 1));
}

bool checkIfMatch(uint32_t leftOli, uint32_t rightOli, uint32_t offset){
    return (leftOli & ((1 << (20-(offset*2))) - 1) ) == ((rightOli >> (offset*2)) & ((1 << (20-(offset*2))) - 1));
}

int main(int argc, char *argv[]) {
    cout << "Test Kondoma" << endl;

    OlinukleoLibrary kondom;
    int startNukleo;
    if (argc > 1) startNukleo = kondom.readFromFile(argv[1],492);
    else startNukleo = kondom.readFromFile("/home/klimas/Documents/Projects/clion/bio/data/negative/113.500-8",492);

    Timer timer;
    timer.start();
    kondom.greedyAlgorithm(startNukleo);
    printf("Processing time: %fs\n", timer.stop());

#ifdef PRINT_STATS
    kondom.printStats();
#endif

    if(!kondom.test())cout << "cos sie zepsulo" << std::endl;

    return 0;
}

void OlinukleoLibrary::greedyAlgorithm(int startNukleo) {
    result.push_back(startNukleo);
    uint32_t successorNucleotide = startNukleo;
    uint32_t predecessorNucleotide = startNukleo;
            cout << "Next " << intIntoStringCoder(startNukleo) << endl;
    int nextNucleotide = -1;
    for (int i = 1; i <= MAX_NEGATIVE; ++i) {
        while( (nextNucleotide = successor(successorNucleotide,i) )!=-1){
//            cout << "Start" <<  intIntoStringCoder(successorNucleotide) << endl;
//            cout << "Next " << intIntoStringCoder(nextNucleotide) << endl;
            successorNucleotide = nextNucleotide;
        };

//    for (int el : result) {
//        cout << intIntoStringCoder(el) << endl;
//    }

        while( (nextNucleotide = predecessor(predecessorNucleotide,i))!=-1){
            predecessorNucleotide = nextNucleotide;
        }
//
//    for (int el : result) {
//        cout << intIntoStringCoder(el) << endl;
//    }

    }

//    for (int el : result) {
//        cout << intIntoStringCoder(el) << endl;
//    }

}

int OlinukleoLibrary::successor(uint32_t nucleotide, unsigned offset) {
    vector<int> successors;
    int base = (nucleotide << (2 * offset)) & ((1 << 20) - 1);
        for (int j = base; j < base + pow(4,offset); ++j) {
            if(microArray[j]) {
#ifdef COUNT_STATS
                ++negativeOffset[offset-1];
#endif
                negativeError+=offset-1;
                successors.push_back(j);
            }
        }
        if(!successors.empty()){
            int r = successors.at(getRandom(0, successors.size() - 1));
            microArray[r] = false;
            result.push_back(r);
            return r;
        }
    return -1;
}

int OlinukleoLibrary::predecessor(uint32_t nucleotide, unsigned offset){
    vector<int> predecessor;
    int step = 262144/pow(4,offset-1);
    int base = (nucleotide>>(2*offset)) & ((1<<20)-1);
        for (int j = base; j <  MATRIX_COUNT; j+=step) {
            if(microArray[j]) {
#ifdef COUNT_STATS
                ++negativeOffset[offset-1];
#endif
                negativeError+=offset-1;
                predecessor.push_back(j);
            }
        }
        if(!predecessor.empty()){
            int r = predecessor.at(getRandom(0, predecessor.size() - 1));
            microArray[r] = false;
            result.push_front(r);
            return r;
        }
    return -1;
}

bool OlinukleoLibrary::test() {
    int i = 0;
    int negative[MAX_NEGATIVE] = {0};
    for (auto element :result){
#ifdef PRINT_RESULT
        cout << intIntoStringCoder(element) << endl;
#endif
        if(i!=result.size()-1) {
            for (int j = 1; j <= MAX_NEGATIVE; ++j) {
                if(checkIfMatch(result.at(i), result.at(i+1),j) ){
                    ++negative[j-1];
                    break;
                }
                if(j==MAX_NEGATIVE) {
                    cout << "nie bylo mnie slychac " << endl;
                            cout << intIntoStringCoder(result.at(i)) << endl;
                    cout << intIntoStringCoder(result.at(i+1)) << endl;
                    return false;
                }
            }
        }
        ++i;
    }
#ifdef PRINT_STATS
    for (int k = 0; k < MAX_NEGATIVE-1; ++k) {
        cout << "Neagtive with " << k+2 << " offset: " << negative[k+1] << endl;
    }
#endif
    return true;
}

void OlinukleoLibrary::printStats() {
    cout << "Słów w sekwencji: " << result.size() << endl;
    cout << "Słów w zbiorze: " << m << endl;
    for (int i = 0; i < MAX_NEGATIVE-1; ++i) {
        cout << "NO " << i + 2 << ": " << negativeOffset[i+1] << endl;
    }
    cout << "Negative: " << negativeError << endl;

}
