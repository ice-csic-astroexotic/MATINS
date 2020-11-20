add_library(grid "")

target_sources(grid
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/grid.f90
    ${CMAKE_CURRENT_LIST_DIR}/ns_structure.f
)

target_link_libraries(grid
  PUBLIC
    constants
    utils
    initial
)