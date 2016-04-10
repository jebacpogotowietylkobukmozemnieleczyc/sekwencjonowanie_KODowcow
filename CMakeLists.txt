cmake_minimum_required(VERSION 2.8.4)
project(bio2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES main.cpp timer.cpp timer.hpp  RandomHeuristic.cpp RandomHeuristic.h Heuristic.cpp Heuristic.h MonteCarlo.cpp MonteCarlo.h)
add_executable(bio2 ${SOURCE_FILES})