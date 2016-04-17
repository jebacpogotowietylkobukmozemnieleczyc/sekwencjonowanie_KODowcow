#include "timer.hpp"
#include "RandomHeuristic.h"
#include "MonteCarlo.h"
#include <iostream>
#include <fstream>
#include <cstring>

void parseAndTest(const char *fname){
    int sub1 = 0, sub2 = 0;
    int oliL = 0, n = 0;
    cout << fname << ";";

    MonteCarlo heuristic;

    //############ NAME PARSER ################
    for( int i = 0; i < strlen(fname); i++){
        if(fname[i] == '.') sub1 = i;
        if(sub1 > 0 && sub2 == 0 && (fname[i] == '-' || fname[i] == '+')) sub2 = i;
    }

    string tmp(fname+sub1+1, fname+sub2);
    oliL = atoi(tmp.c_str());

    n = oliL;
    cout << n << ";";

    string tmp2(fname+sub2+1, fname+strlen(fname));
    if(fname[sub2] == '-') oliL -= atoi(tmp2.c_str());
    else oliL += atoi(tmp2.c_str());
    //########### NAME PARSER END ###########

    heuristic.readFromFile(fname, oliL, n);

    heuristic.initRandomVector();

    Timer timer;
    timer.start();

    heuristic.run();

    printf("%f;", timer.stop());

    heuristic.printStats();

    cout << heuristic.checkIfCorrect() << endl; //sprawdza czy się mieści

    tmp.clear();
    tmp2.clear();

}

void parseAndTestAll(const char *fname){
    fstream file;
    int iterator = 0;
    string tmp;

    file.open(fname, ios::in);
    while (!file.eof()){
        file >> tmp;
        const char * c = tmp.c_str();
        parseAndTest(c);
        iterator++;
    }
    file.close();
}


int main(int argc, char *argv[]) {
    cout << "Test" << endl;

    parseAndTestAll("/home/klimas/Documents/Projects/clion/bio2/testInstances.txt");

}

