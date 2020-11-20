add_library(initial "")

target_sources(initial
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/input.f90
    ${CMAKE_CURRENT_LIST_DIR}/binit.f90
    ${CMAKE_CURRENT_LIST_DIR}/resume.f90
)

target_link_libraries(initial
  PUBLIC
    bevol
    grid
    constants
    math
)

