cmake_minimum_required(VERSION 3.18)
project(pyi2em LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED gsl)

add_subdirectory(src)