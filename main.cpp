#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <random>
#include <bits/stl_deque.h>

#define A_BIN 0
#define C_BIN 1
#define G_BIN 2
#define T_BIN 3

#define OLI_LENGTH 10
#define MATRIX_COUNT 1048576


#define MAX_NEGATIVE 3
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
    deque<int> result;
    unsigned n;

public:
    bool microArray[MATRIX_COUNT] = {false};
    uint32_t stringIntoIntCoder(string sequence);
    string intIntoStringCoder(uint32_t codedNumber);
    int readFromFile(const char* fname, unsigned fileLength);
    void printVector();
    void printVectorAsString();

    void greedyAlgorithm(int startNukleo);
    void successor(int Nukleo);
    void predecessor(int Nukleo);


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
//        oliLib.push_back(OlinukleoLibrary::stringIntoIntCoder(tmp));
        microArray[OlinukleoLibrary::stringIntoIntCoder(tmp)]=true;
        //todo differFileLength
        if(iterator==randomLine)randomNukleo=OlinukleoLibrary::stringIntoIntCoder(tmp);
        iterator++;
    }

    file.close();

    n = iterator;
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
    if ( (leftOli & ((1<<18)-1) ) == ( (rightOli>>2) & ((1<<18)-1)) ) return true; else return false;
}


int main(int argc, char *argv[]) {
    cout << "Test Kondoma" << endl;

    OlinukleoLibrary kondom;
    int startNukleo;
    if (argc > 1) startNukleo = kondom.readFromFile(argv[1],500);
    else startNukleo = kondom.readFromFile("/home/klimas/Documents/Projects/clion/bio/data/negative/113.500-8",500);
    kondom.greedyAlgorithm(startNukleo);


    //kondom.printVectorAsString();

//    cout << getRandom(0,7) << endl;
//    cout << kondom.microArray[549385] << endl;

    return 0;
}

void OlinukleoLibrary::greedyAlgorithm(int startNukleo) {
    result.push_back(startNukleo);
//    successor(startNukleo);
    predecessor(startNukleo);

    for (auto element :result){
        cout << intIntoStringCoder(element) << endl;
    }


    cout << checkIfMatch(result.at(0), result.at(1)) << endl;
}

void OlinukleoLibrary::successor(int Nukleo) {
    int base;
    for (int i = 1; i <= MAX_NEGATIVE; ++i) {
        base = (Nukleo<<2) & ((1<<20)-1);
        for (int j = base; j < base + pow(4,MAX_NEGATIVE); ++j) {
            if(microArray[j]) {
                result.push_back(j);
//                successor(j);
                return;

            }
        }

    }

}

void OlinukleoLibrary::predecessor(int Nukleo){
    int base;
    int step = 262144;
    for (int i = 1; i <= MAX_NEGATIVE; ++i) {
        base = (Nukleo>>2) & ((1<<20)-1);
        for (int j = base; j < base + MATRIX_COUNT; j+=step) {
            if(microArray[j]) {
                result.push_front(j);
//                predecessor(j);
                return;
            }
        }
        step/=4;

    }

}
