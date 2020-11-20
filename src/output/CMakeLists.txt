add_library(output "")

target_link_libraries(output
  PUBLIC
    constants
    utils
    bevol
    grid
)

target_sources(output
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/output.f90
)
