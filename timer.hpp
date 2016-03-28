#ifndef TIMER_HPP
#define TIMER_HPP

#ifdef __WIN32
    #include <windows.h>
#else
    #include <time.h>
#endif // __WIN32

class Timer {
public:
    Timer();

    void start();
    double stop();

private:
    double lastTime;

#ifdef __WIN32
    LARGE_INTEGER freq;
    LARGE_INTEGER t0;
#else
    struct timespec t0;
#endif // __WIN32
};

#endif // TIMER_HPP
