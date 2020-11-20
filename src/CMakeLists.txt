add_executable(mt2d
  main.f90
)

# Create a list of internal library link targets to link with MT2D.
# Make sure to add each one of the required components to MT2D_LIBS.
# Note: AMBIPOLAR and MAGNETOSPHERE can be disabled.

option(ENABLE_AMBIPOLAR "Enable Ambipolar" OFF)
option(ENABLE_MAGNETOSPHERE "Enable Magnetosphere" OFF)

list(APPEND MT2D_LIBS "")

if(ENABLE_AMBIPOLAR)
  add_subdirectory(ambipolar)
  list(APPEND MT2D_LIBS "ambipolar")
  message("Compiling with AMBIPOLAR as requested")
endif()

add_subdirectory(bevol)
list(APPEND MT2D_LIBS "bevol")

add_subdirectory(constants)
list(APPEND MT2D_LIBS "constants")

add_subdirectory(thermal_evol)
list(APPEND MT2D_LIBS "thermal_evol")

add_subdirectory(math)
list(APPEND MT2D_LIBS "math")

if(ENABLE_MAGNETOSPHERE)
  add_subdirectory(magnetosphere)
  list(APPEND MT2D_LIBS "magnetosphere")
  message("Compiling with MAGNETOSPHERE as requested")
endif()

add_subdirectory(initial)
list(APPEND MT2D_LIBS "initial")

add_subdirectory(grid)
list(APPEND MT2D_LIBS "grid")

add_subdirectory(microphys)
list(APPEND MT2D_LIBS "microphys")

add_subdirectory(output)
list(APPEND MT2D_LIBS "output")

add_subdirectory(utils)
list(APPEND MT2D_LIBS "utils")

target_link_libraries(mt2d
  PRIVATE
    ${MT2D_LIBS}
)

# Custom commands to execute before of after compilation.
#   * Create build/bin/in directory.
#   * Create build/bin/out directory.
#   * Create build/bin/outb directory.
#   * Copy in/* to build/bin/in

add_custom_command(
        TARGET mt2d PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_BINARY_DIR}/bin/in)
add_custom_command(
        TARGET mt2d PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_BINARY_DIR}/bin/out)
add_custom_command(
        TARGET mt2d PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_BINARY_DIR}/bin/outb)

file(COPY ${CMAKE_SOURCE_DIR}/in/
     DESTINATION ${CMAKE_BINARY_DIR}/bin/in/)