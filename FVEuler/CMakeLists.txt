cmake_minimum_required(VERSION 3.25)
project(FVEuler)

set(CMAKE_CXX_STANDARD 17)

add_executable(FVEuler main.cpp EulerFlux.cpp EulerFlux.h BoundaryConditions.cpp BoundaryConditions.h FileIO.cpp FileIO.h Indexing.h)
