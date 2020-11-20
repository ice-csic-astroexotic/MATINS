add_library(math "")

target_link_libraries(math
  PUBLIC
    constants
    grid
)

target_sources(math
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/math.f90
    ${CMAKE_CURRENT_LIST_DIR}/legpol.f90
    ${CMAKE_CURRENT_LIST_DIR}/odeint.f
)
