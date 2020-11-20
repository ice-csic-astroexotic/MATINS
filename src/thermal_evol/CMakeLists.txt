add_library(thermal_evol "")

target_sources(thermal_evol
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/fluxes_boundary.f90
    ${CMAKE_CURRENT_LIST_DIR}/tevol.f90
    ${CMAKE_CURRENT_LIST_DIR}/timestep.f90
)

target_link_libraries(thermal_evol
  PUBLIC
    bevol
    constants
    microphys
    math
)