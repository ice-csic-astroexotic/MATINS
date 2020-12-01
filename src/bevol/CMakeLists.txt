add_library(bevol "")

target_sources(bevol
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/magnetic_evol.f90
    ${CMAKE_CURRENT_LIST_DIR}/break.f90
    ${CMAKE_CURRENT_LIST_DIR}/magnetic_analysis.f90
    ${CMAKE_CURRENT_LIST_DIR}/pevol.f90
    ${CMAKE_CURRENT_LIST_DIR}/btimestep.f90
)

target_link_libraries(bevol
  PUBLIC
    constants
    grid
    initial
    math
    output
)

