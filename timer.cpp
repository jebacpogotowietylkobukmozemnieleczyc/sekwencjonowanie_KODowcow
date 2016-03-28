#include "timer.hpp"

#ifdef __WIN32

Timer::Timer() {
    lastTime = 0;
    QueryPerformanceFrequency(&freq);
}

void Timer::start() {
    QueryPerformanceCounter(&t0);
}

double Timer::stop() {
    LARGE_INTEGER t1;
    QueryPerformanceCounter(&t1);

    lastTime = ((double)(t1.QuadPart - t0.QuadPart)) / freq.QuadPart;

    return lastTime;
}

#else

Timer::Timer() {
    lastTime = 0;
}

void Timer::start() {
    clock_gettime(CLOCK_MONOTONIC, &t0);
}

double Timer::stop() {
    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);

    lastTime = t1.tv_sec - t0.tv_sec + ((double)(t1.tv_nsec - t0.tv_nsec))
               / 1000000000;

    return lastTime;
}

#endif // __WIN32
