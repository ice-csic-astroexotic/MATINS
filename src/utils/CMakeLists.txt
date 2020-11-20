add_library(utils "")

target_link_libraries(utils
  PUBLIC
)

target_sources(utils
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/utils.f90
)
