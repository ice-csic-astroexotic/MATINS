cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(MAGNETOTHERMAL2D LANGUAGES Fortran C)

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Set custom flags for the different build types.
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -ffpe-trap=invalid,zero -ffixed-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,zero -fcheck=all -Wall -ffixed-line-length-none -pg")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -ffpe-trap=invalid,zero -ffixed-line-length-none -pg")

# Setting default build type to Release.
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# The main subdirectory that contains the target and libraries.
add_subdirectory(src)

# Contains an example "external" library we will link to.
add_subdirectory(external)

# Add tests.
if (ENABLE_TEST)
  message(STATUS "Tests enabled! So building tests for ya!")
  enable_testing()
  add_subdirectory(tests)
endif()

# Add benchmarks
if (ENABLE_BENCHMARKS)
  message(STATUS "Benchmarks enabled! Let me compile them for ya!")
  add_subdirectory(benchmarks)
endif()
