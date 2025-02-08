cmake_minimum_required(VERSION 3.29)
project(mcccm_fortran LANGUAGES Fortran C CXX)

set(AUX_SOURCES
        generated_geometry.h
        geometry.h
        generated_geometry.cpp
        geometry.cpp
        geometry_c_binding.cpp
        geometry_c_binding.h
        geometry.f90
)

add_executable(mcccm_fortran
        main.f90
        ${AUX_SOURCES}
)

add_executable(mcccm_fortran_tests
        test/unit_tests.f90
        ${AUX_SOURCES}
)
