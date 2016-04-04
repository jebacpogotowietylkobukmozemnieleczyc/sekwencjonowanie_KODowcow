//
// Created by klimas on 04.04.16.
//

#include "RandomHeuristic.h"


void RandomHeuristic::runRandom() {
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

bool RandomHeuristic::test() {
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

void RandomHeuristic::printStats() {
    cout << "Słów w sekwencji: " << result.size() << endl;
    cout << "Słów w zbiorze: " << m << endl;
    for (int i = 0; i < MAX_NEGATIVE-1; ++i) {
        cout << "NO " << i + 2 << ": " << countOffset[i+1] << endl;
    }
    cout << "Negative: " << negativeError << endl;
}

void RandomHeuristic::shuffleVector(int offset) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(randomVectors[offset - 1].begin(), randomVectors[offset - 1].end(), g);
}

void RandomHeuristic::generateRandomVector(unsigned min, unsigned max) {
    unsigned s = (min + max) / 2 * (max - min + 1);
    unsigned n = m/s;
    int value(0);
    std::generate_n(std::back_inserter(randomVector), n+1, [value]()mutable { return value++; });

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(randomVector.begin(), randomVector.end(), g);
}

void RandomHeuristic::initRandomVector() {


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

void RandomHeuristic::run() {
    runRandom();
}


uint32_t RandomHeuristic::stringIntoIntCoder(string sequence){
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


string RandomHeuristic::intIntoStringCoder(uint32_t codedNumber){
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

void RandomHeuristic::readFromFile(const char *fname, unsigned fileLength){
    fstream file;
    int iterator = 0;
    string tmp;
    int randomNucleotide = -1;

    int randomLine = getRandom(0, fileLength - 1);
    file.open(fname, ios::in);
    while (!file.eof()){
        file >> tmp;
        microArray[RandomHeuristic::stringIntoIntCoder(tmp)]=true;
        //todo take properties from file name
        if(iterator==randomLine){
            randomNucleotide = RandomHeuristic::stringIntoIntCoder(tmp);
            microArray[RandomHeuristic::stringIntoIntCoder(tmp)]= false;
        }
        iterator++;
    }

    file.close();

    //todo empty line at end
    m = iterator - 1;
    result.push_back(randomNucleotide);
}

void RandomHeuristic::printResult(){
    for (int element : result) {
        cout << element << endl;

    }
}

void RandomHeuristic::printResultAsString(){
    for (int element :result){
        cout << RandomHeuristic::intIntoStringCoder(element) << endl;
    }
}



bool RandomHeuristic::checkIfMatch(uint32_t leftNucleotide, uint32_t rightNucleotide) {
    return (leftNucleotide & ((1 << 18) - 1) ) == ((rightNucleotide >> 2) & ((1 << 18) - 1));
}


bool RandomHeuristic::checkIfMatch(uint32_t leftNucleotide, uint32_t rightNucleotide, uint32_t offset) {
    return (leftNucleotide & ((1 << (20 - (offset * 2))) - 1) ) == ((rightNucleotide >> (offset * 2)) & ((1 << (20 - (offset * 2))) - 1));
}

uint32_t RandomHeuristic::getRandom(uint32_t min, uint32_t max){
std::uniform_int_distribution<mt19937::result_type > udist(min, max);
mt19937 rand;
random_device randomDevice;
mt19937::result_type const seedval = randomDevice(); // get this from somewhere
rand.seed(seedval);

mt19937::result_type random_number = udist(rand);

return random_number;
}
