# bindings.cpp should create the py module and include core sources
include(FetchContent)

# Fetch cubature at build time (no system dep required for wheels)
# Use git tag to avoid tarball hash drift on GitHub auto-generated archives.
# Force a static build so our Python extension has no runtime dependency on
# libcubature.so/.dylib (fixes ImportError in tests and simplifies wheels).
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
FetchContent_Declare(
  cubature
  GIT_REPOSITORY https://github.com/stevengj/cubature.git
  GIT_TAG v1.0.4
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cubature)

# Ensure the static library is built with -fPIC so it can be linked into
# our Python extension module on Linux.
if(TARGET cubature)
  set_target_properties(cubature PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

pybind11_add_module(pyi2em
  bindings.cpp
  em_i2em.cpp)

target_include_directories(pyi2em PRIVATE
  ${GSL_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/src
  ${cubature_SOURCE_DIR})

target_link_directories(pyi2em PRIVATE ${GSL_LIBRARY_DIRS})
target_link_libraries(pyi2em PRIVATE ${GSL_LIBRARIES} cubature)

if(NOT MSVC)
  target_compile_options(pyi2em PRIVATE -fvisibility=hidden)
endif()

# Install the extension at the top-level of the wheel (importable as `pyi2em`)
install(TARGETS pyi2em
  LIBRARY DESTINATION .
  ARCHIVE DESTINATION .
  RUNTIME DESTINATION .)
