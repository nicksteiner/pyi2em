cmake_minimum_required(VERSION 3.0.0)
project(em_i2em VERSION 0.1.0)

include(CTest)
enable_testing()

set (CMAKE_CXX_STANDARD 11)

find_library(GSL_LIBRARY      NAMES gsl)
find_library(GSL_LIBRARY_BLAS NAMES gslcblas)

include_directories(lib/cubature-1.0 lib/gsl src)

set(PYBIND11_CPP_STANDARD -std=c++11)
add_subdirectory(pybind11)

set(SOURCES src/main.cpp lib/cubature-1.0/hcubature.c lib/cubature-1.0/pcubature.c src/em_i2em.cpp)

pybind11_add_module(pyi2em SHARED ${SOURCES})

target_link_libraries(pyi2em PRIVATE ${GSL_LIBRARY} ${GSL_LIBRARY_BLAS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)