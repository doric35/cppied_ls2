cmake_minimum_required(VERSION 3.5)
project(algorithms)
message(STATUS "CMAKE SOUCRE DIR: " ${CMAKE_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(local_search_algorithms
        main/local_search_main.cpp
        local_search_algorithms/ls2.cpp
        cppied_algorithm.cpp)

add_executable(initialization_algorithms
        main/initialization_main.cpp
        initialization_algorithms/dp_sweeper.cpp
        cppied_algorithm.cpp)
