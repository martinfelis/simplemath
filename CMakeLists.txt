cmake_minimum_required(VERSION 3.12)
project(simplemathv2)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(include/SimpleMath)
include_directories(tests)

add_executable(simplemathv2
        include/SimpleMath/SimpleMath.h
        tests/main.cpp
        tests/tests.cpp)
