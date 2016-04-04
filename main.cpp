#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <random>
#include <array>
#include <algorithm>

#include "timer.hpp"

#define A_BIN 0
#define C_BIN 1
#define G_BIN 2
#define T_BIN 3

#define OLI_LENGTH 10
#define MATRIX_COUNT 1048576
#define MAX_NEGATIVE 9

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
    OlinukleoLibrary() {


    }

private:
    array< vector<int>,MAX_NEGATIVE> randomVectors;
    vector<unsigned > randomVector;
    unsigned randomVectorIterator;
    array<unsigned,MAX_NEGATIVE> randomVectorLimit;

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

    void greedyAlgorithm(int startNukleo, unsigned randomFrequency);
    int successor(uint32_t base, unsigned offset);
    int predecessor(uint32_t nucleotide, unsigned offset);
    void shuffleVector(int offset);
    void generateRandomVector(unsigned min, unsigned max);
    void init();


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

    kondom.init();

    Timer timer;
    timer.start();
    kondom.greedyAlgorithm(startNukleo,10);
    printf("Processing time: %fs\n", timer.stop());

#ifdef PRINT_STATS
    kondom.printStats();
#endif

    if(!kondom.test())cout << "cos sie zepsulo" << std::endl;

    return 0;
}

void OlinukleoLibrary::greedyAlgorithm(int startNukleo,unsigned randomFrequency) {
    result.push_back(startNukleo);
    uint32_t successorNucleotide = startNukleo;
    uint32_t predecessorNucleotide = startNukleo;

    array<int,MAX_NEGATIVE> countOffset = {0};


    int nextNucleotide = -1;
    for (int i = 1; i <= MAX_NEGATIVE; ++i) {
        while( (nextNucleotide = successor(successorNucleotide,i) )!=-1){
//            cout << "Start" <<  intIntoStringCoder(successorNucleotide) << endl;
//            cout << "Next " << intIntoStringCoder(nextNucleotide) << endl;
            if( (++countOffset[i-1])==randomVectorLimit[i-1]){
                shuffleVector(i);
                randomVectorIterator = randomVectorIterator >= randomVector.size() ? 0 : randomVectorIterator;
                randomVectorLimit[i-1]+=randomVector[randomVectorIterator++];
            }
            if(i>1)--i;
            successorNucleotide = nextNucleotide;
        };

//    for (int el : result) {
//        cout << intIntoStringCoder(el) << endl;
//    }

        while( (nextNucleotide = predecessor(predecessorNucleotide,i))!=-1){
            if( (++countOffset[i-1])==randomVectorLimit[i-1]){
                shuffleVector(i);
                randomVectorIterator = randomVectorIterator >= randomVector.size() ? 0 : randomVectorIterator;
                randomVectorLimit[i-1]+=randomVector[randomVectorIterator++];
            }
            if(i>1)--i;
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
    int base = (nucleotide << (2 * offset)) & ((1 << 20) - 1);
        for (int j = 0; j <  pow(4,offset); ++j) {
            int r = base + randomVectors[offset - 1].at(j);
            if(microArray[r]) {
#ifdef COUNT_STATS
                ++negativeOffset[offset-1];
#endif
                negativeError+=offset-1;

                result.push_back(r);
                microArray[r] = false;
                return r;
            }
        }
    return -1;
}

int OlinukleoLibrary::predecessor(uint32_t nucleotide, unsigned offset){
    int step = 262144/pow(4,offset-1);
    int base = (nucleotide>>(2*offset)) & ((1<<20)-1);
    for (int j = 0; j <  pow(4,offset); ++j) {
        int r = base + randomVectors[offset - 1].at(j) * step;
            if(microArray[j]) {
#ifdef COUNT_STATS
                ++negativeOffset[offset-1];
#endif
                negativeError+=offset-1;

                microArray[r] = false;
                result.push_front(r);
                return r;
            }
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

void OlinukleoLibrary::shuffleVector(int offset) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(randomVectors[offset - 1].begin(), randomVectors[offset - 1].end(), g);
}

void OlinukleoLibrary::generateRandomVector(unsigned min, unsigned max) {
    unsigned s = (min + max) / 2 * (max - min + 1);
    unsigned n = m/s;
    int value(0);
    std::generate_n(std::back_inserter(randomVector), n+1, [value]()mutable { return value++; });

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(randomVector.begin(), randomVector.end(), g);
}

void OlinukleoLibrary::init() {


    //wygenerowanie wektorów przejsc
    for (int i = 0; i < MAX_NEGATIVE; ++i) {
        std::generate_n(std::back_inserter(randomVectors[i]), pow(4, i + 1), [&](){ return randomVectors[i].size(); });
        shuffleVector(i + 1);
//            cout << "Array: " <<  i << endl;
//            for (int el : randomVectors[i]) {
//                cout << el << endl;
//            }
    }

    //wygenerowanie wektorow okreslajacych co ile wektory przejsc beda "szaflowane"
    generateRandomVector(20, 40);

    //przypisanie poczatkowych wartości do limitu przy ktorym vectory przejsc beda "szaflowane"
    for (int j = 0; j < MAX_NEGATIVE; ++j) {
        randomVectorLimit[j] = randomVectors[j].at(0);
    }


}
