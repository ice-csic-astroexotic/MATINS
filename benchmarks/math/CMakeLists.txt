add_executable(benchmark_math_ludcmp
  benchmark_math_ludcmp.f90
)

target_link_libraries(benchmark_math_ludcmp
  PUBLIC
    constants
    math
)

set_target_properties(benchmark_math_ludcmp
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks"
)

###############################################################

add_executable(benchmark_math_getugam
  benchmark_math_getugam.f90
)

target_link_libraries(benchmark_math_getugam
  PUBLIC
    constants
    math
)

set_target_properties(benchmark_math_getugam
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks"
)

###############################################################

add_executable(benchmark_math_compute_gausslegendre_weights
  benchmark_math_compute_gausslegendre_weights.f90
)

target_link_libraries(benchmark_math_compute_gausslegendre_weights
  PUBLIC
    constants
    math
)

set_target_properties(benchmark_math_compute_gausslegendre_weights
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks"
)