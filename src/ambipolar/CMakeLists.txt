add_library(ambipolar "")

target_sources(ambipolar
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/ambipolar_interface.f
    ${CMAKE_CURRENT_LIST_DIR}/ambipolar_sf.f
    ${CMAKE_CURRENT_LIST_DIR}/ambipolar.f
    ${CMAKE_CURRENT_LIST_DIR}/trdig.f
)

target_link_libraries(ambipolar
  PUBLIC
    bevol
    grid
    math)
