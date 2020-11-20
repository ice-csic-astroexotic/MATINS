add_library(microphys "")

target_link_libraries(microphys
  PUBLIC
    constants
    grid
    initial
)

target_sources(microphys
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/microphys_interface.f90
    ${CMAKE_CURRENT_LIST_DIR}/cond.f
    ${CMAKE_CURRENT_LIST_DIR}/condegin.f
    ${CMAKE_CURRENT_LIST_DIR}/conduct19.f
    ${CMAKE_CURRENT_LIST_DIR}/eff_mass.f
    ${CMAKE_CURRENT_LIST_DIR}/emissivity.f
    ${CMAKE_CURRENT_LIST_DIR}/eos14.f
    ${CMAKE_CURRENT_LIST_DIR}/eosmag14.f
    ${CMAKE_CURRENT_LIST_DIR}/potekhin2019.f
    ${CMAKE_CURRENT_LIST_DIR}/superfluid.f
    ${CMAKE_CURRENT_LIST_DIR}/cvf.f90
    ${CMAKE_CURRENT_LIST_DIR}/shear.f90
    ${CMAKE_CURRENT_LIST_DIR}/get_Zimp.f90
    ${CMAKE_CURRENT_LIST_DIR}/cphonon_2010.f
    ${CMAKE_CURRENT_LIST_DIR}/cphonon_Chu.f
)
