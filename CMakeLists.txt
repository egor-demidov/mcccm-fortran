cmake_minimum_required(VERSION 3.29)
project(mcccm_fortran LANGUAGES Fortran C CXX)

add_executable(mcccm_fortran
        main.f90
        generated_geometry.h
        geometry.h
        generated_geometry.cpp
        geometry.cpp
        geometry_c_binding.cpp
        geometry_c_binding.h
        geometry.f90
)
