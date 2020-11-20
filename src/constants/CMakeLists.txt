add_library(constants "")

target_sources(constants
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/reals.f90
    ${CMAKE_CURRENT_LIST_DIR}/constants.f90
)
