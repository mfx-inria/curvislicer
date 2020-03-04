cmake_minimum_required(VERSION 2.6)
project(curvislice)

add_subdirectory(libs/SolverWrapper)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}/libs/LibSL-small/src
        ${PROJECT_SOURCE_DIR}/libs/LibSL-small/src/LibSL
        ${PROJECT_SOURCE_DIR}/libs/tclap/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/src/
)

CONFIGURE_FILE(
        "${CMAKE_CURRENT_SOURCE_DIR}/src/path.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
)

SET(LIBSL_SOURCES
        libs/LibSL-small/src/LibSL/Math/Math.cpp
        libs/LibSL-small/src/LibSL/Math/Vertex.cpp
        libs/LibSL-small/src/LibSL/Mesh/Mesh.cpp
        libs/LibSL-small/src/LibSL/Mesh/MeshFormat_stl.cpp
        libs/LibSL-small/src/LibSL/Mesh/VertexFormat_dynamic.cpp
        libs/LibSL-small/src/LibSL/System/System.cpp
        libs/LibSL-small/src/LibSL/Image/Image.cpp
        libs/LibSL-small/src/LibSL/CppHelpers/CppHelpers.cpp
        libs/LibSL-small/src/LibSL/StlHelpers/StlHelpers.cpp
        libs/LibSL-small/src/loki/NullType.h
        libs/LibSL-small/src/loki/static_check.h
        libs/LibSL-small/src/loki/Typelist.h
        libs/LibSL-small/src/loki/TypelistMacros.h
        libs/LibSL-small/src/loki/TypeManip.h
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
        ${LIBSL_SOURCES}
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

target_link_libraries(curvislice_osqp SolverWrapper)

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
    ${LIBSL_SOURCES}
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
  
  target_link_libraries(curvislice_grb SolverWrapper)
  
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
  ${LIBSL_SOURCES}
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

target_link_libraries(uncurve SolverWrapper)

if(BUILD_WITH_GRB)
  target_compile_options(uncurve PRIVATE -DGRB)
endif(BUILD_WITH_GRB)

####################################

ADD_EXECUTABLE(reorient
  src/reorient.cpp
  ${LIBSL_SOURCES}
)
INSTALL(
  TARGETS reorient
  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)
set_target_properties(reorient
  PROPERTIES
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF
)

####################################

ADD_EXECUTABLE(subdiv
  src/subdiv.cpp
  ${LIBSL_SOURCES}
)
INSTALL(
  TARGETS subdiv
  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)

####################################


if(WIN32)
  if (BUILD_WITH_GRB)
    target_link_libraries(curvislice_grb shlwapi)
  endif(BUILD_WITH_GRB)
  target_link_libraries(curvislice_osqp shlwapi)
  target_link_libraries(reorient shlwapi)
  target_link_libraries(subdiv shlwapi)
  target_link_libraries(uncurve shlwapi)
endif(WIN32)