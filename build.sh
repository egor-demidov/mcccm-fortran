#!/bin/bash

# Create a build directory for an out-of-tree build
mkdir -p build
cd build

# Compile the Fortran sources
gfortran -c ../geometry.f90
gfortran -c ../main.f90
# Compile the C++ sources
g++ -c ../geometry.cpp -o geometry_cxx.o    # custom object file name to prevent name overlaps
g++ -c ../generated_geometry.cpp
g++ -c ../geometry_c_binding.cpp
# Link the object files into an executable (-lgfortran is required to link the Fortran standard library)
g++ *.o -lgfortran -o main

# Return to source directory
cd ..
