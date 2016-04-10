
#include "timer.hpp"
#include "RandomHeuristic.h"
#include "MonteCarlo.h"


int main(int argc, char *argv[]) {
    cout << "Test" << endl;

    MonteCarlo heuristic;
    if (argc > 1)  heuristic.readFromFile(argv[1], 492);
    else  heuristic.readFromFile("/home/klimas/Documents/Projects/clion/bio/data/negative/113.500-8", 492);

    heuristic.initRandomVector();

    Timer timer;
    timer.start();
    heuristic.run();
    printf("Processing time: %fs\n", timer.stop());

#ifdef PRINT_STATS
    heuristic.printStats();
#endif

    if(!heuristic.test())cout << "cos sie zepsulo" << std::endl;

    return 0;
}

