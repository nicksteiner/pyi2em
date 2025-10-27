# bindings.cpp should create the py module and include core sources
pybind11_add_module(pyi2em
  bindings.cpp
  em_i2em.cpp
  ../lib/cubature-1.0/hcubature.c
  ../lib/cubature-1.0/pcubature.c
  ../lib/cubature-1.0/clencurt_gen.c)

target_include_directories(pyi2em PRIVATE
  ${GSL_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/lib/cubature-1.0)

target_link_directories(pyi2em PRIVATE ${GSL_LIBRARY_DIRS})
target_link_libraries(pyi2em PRIVATE ${GSL_LIBRARIES})

if(NOT MSVC)
  target_compile_options(pyi2em PRIVATE -fvisibility=hidden)
endif()
