cmake_minimum_required(VERSION 3.14)
project(potts)

set(CMAKE_CXX_STANDARD 14)

add_executable(potts main.cpp some_functions.cpp)
target_link_libraries (potts armadillo)
