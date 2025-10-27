# bindings.cpp should create the py module and include core sources
include(FetchContent)

# Fetch cubature at build time (no system dep required for wheels)
FetchContent_Declare(
  cubature
  URL https://github.com/stevengj/cubature/archive/refs/tags/v1.0.4.tar.gz
  URL_HASH SHA256=8f3fc3b2fa5474f9c6d3f2c8374e036f7a205b786c396f26faa9c67806bfe5a6
)
FetchContent_MakeAvailable(cubature)

set(CUBATURE_SRC
  ${cubature_SOURCE_DIR}/hcubature.c
  ${cubature_SOURCE_DIR}/pcubature.c
  ${cubature_SOURCE_DIR}/clencurt_gen.c)

pybind11_add_module(pyi2em
  bindings.cpp
  em_i2em.cpp
  ${CUBATURE_SRC})

target_include_directories(pyi2em PRIVATE
  ${GSL_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/src
  ${cubature_SOURCE_DIR})

target_link_directories(pyi2em PRIVATE ${GSL_LIBRARY_DIRS})
target_link_libraries(pyi2em PRIVATE ${GSL_LIBRARIES})

if(NOT MSVC)
  target_compile_options(pyi2em PRIVATE -fvisibility=hidden)
endif()
