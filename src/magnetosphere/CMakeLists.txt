add_library(magnetosphere "")

target_sources(magnetosphere
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/energy.f
    ${CMAKE_CURRENT_LIST_DIR}/force_free.f
    ${CMAKE_CURRENT_LIST_DIR}/interface.f
    ${CMAKE_CURRENT_LIST_DIR}/multipole.f
    ${CMAKE_CURRENT_LIST_DIR}/operator.f
    ${CMAKE_CURRENT_LIST_DIR}/toroidal_fn.f
    ${CMAKE_CURRENT_LIST_DIR}/twist.f
    ${CMAKE_CURRENT_LIST_DIR}/surface_transition.f
    ${CMAKE_CURRENT_LIST_DIR}/covsrt.f
    ${CMAKE_CURRENT_LIST_DIR}/gaussj.f
    ${CMAKE_CURRENT_LIST_DIR}/legendre.f
    ${CMAKE_CURRENT_LIST_DIR}/lfit.f
    ${CMAKE_CURRENT_LIST_DIR}/mrqmin.f
    ${CMAKE_CURRENT_LIST_DIR}/mrqmin2d.f
)
