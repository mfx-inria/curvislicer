####################################

cmake_minimum_required(VERSION 2.6)
project(curvislice)

add_subdirectory(libs/SolverWrapper)
add_subdirectory(libs/LibSL-small)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}/libs/tclap/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/src/
)

####################################

ADD_EXECUTABLE(curvislice_osqp
        src/main.cpp
        src/gcode.h
        src/gcode.cpp
        src/helpers.h
        src/TetMesh.h
        src/TetMesh.cpp
        src/thicknesses.h
        src/MeshFormat_msh.h
        src/MeshFormat_msh.cpp
        )

INSTALL(
  TARGETS curvislice_osqp
  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)

set_target_properties(curvislice_osqp
  PROPERTIES
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF
)

target_link_libraries(curvislice_osqp SolverWrapper LibSL)

target_compile_options(curvislice_osqp PRIVATE -DOSQP)

####################################

if (BUILD_WITH_GRB)
  ADD_EXECUTABLE(curvislice_grb
    src/main.cpp
    src/gcode.h
    src/gcode.cpp
    src/helpers.h
    src/TetMesh.h
    src/TetMesh.cpp
    src/thicknesses.h
    src/MeshFormat_msh.h
    src/MeshFormat_msh.cpp
  )

  INSTALL(
    TARGETS curvislice_grb
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
  )

  set_target_properties(curvislice_grb
    PROPERTIES
      CXX_STANDARD 17
      CXX_EXTENSIONS OFF
  )

  target_link_libraries(curvislice_grb SolverWrapper LibSL)

  target_compile_options(curvislice_grb PRIVATE -DGRB)
  target_compile_definitions(curvislice_grb PUBLIC "-DHAS_GUROBI")
endif(BUILD_WITH_GRB)

####################################

ADD_EXECUTABLE(uncurve
  src/uncurve.cpp
  src/gcode.h
  src/gcode.cpp
  src/TetMesh.h
  src/TetMesh.cpp
  src/thicknesses.h
  src/MeshFormat_msh.h
  src/MeshFormat_msh.cpp
)

INSTALL(
  TARGETS uncurve
  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)

set_target_properties(uncurve
  PROPERTIES
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF
)

target_link_libraries(uncurve SolverWrapper LibSL)

if(BUILD_WITH_GRB)
  target_compile_options(uncurve PRIVATE -DGRB)
else(BUILD_WITH_GRB)
  target_compile_options(uncurve PRIVATE -DOSQP)
endif(BUILD_WITH_GRB)

####################################

if(WIN32)
  if (BUILD_WITH_GRB)
    target_link_libraries(curvislice_grb shlwapi)
  endif(BUILD_WITH_GRB)
  target_link_libraries(curvislice_osqp shlwapi)
  target_link_libraries(uncurve shlwapi)
endif(WIN32)
