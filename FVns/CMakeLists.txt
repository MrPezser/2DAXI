cmake_minimum_required(VERSION 3.25)
project(FVns)

set(CMAKE_CXX_STANDARD 17)

add_executable(FVns main.cpp EulerFlux.cpp EulerFlux.h BoundaryConditions.cpp BoundaryConditions.h FileIO.cpp FileIO.h Indexing.h SpatialDiscretization.cpp SpatialDiscretization.h MeshModule.cpp MeshModule.h StateVariables.h LUtools.cpp LUtools.h Jacobian.cpp Jacobian.h)
