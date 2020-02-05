cmake_minimum_required(VERSION 2.6)
project(curvislice)

#OPTION(CURVI_BUILD_WITH_OSPQ  "Build with OSQ (experimental, work in progress)" OFF)

#if (CURVI_BUILD_WITH_OSPQ)
#add_subdirectory(libs/osqp)
#add_definitions(-DUSE_OSQP)
#set(BUILD_WITH_OSPQ 1)
#endif()

#add_subdirectory(libs/SolverWrapper)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}/LibSL-small/src
        ${PROJECT_SOURCE_DIR}/LibSL-small/src/LibSL
        $ENV{GUROBI_HOME}/include
        ${PROJECT_SOURCE_DIR}/libs/tclap/include/
#        ${PROJECT_SOURCE_DIR}/libs/osqp/include/
#        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/include/
#        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/src/
)

link_directories(
        $ENV{GUROBI_HOME}/lib/
)

CONFIGURE_FILE(
        "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
)

SET(LIBSL_SOURCES
        LibSL-small/src/LibSL/Math/Math.cpp
        LibSL-small/src/LibSL/Math/Vertex.cpp
        LibSL-small/src/LibSL/Mesh/Mesh.cpp
        LibSL-small/src/LibSL/Mesh/MeshFormat_stl.cpp
        LibSL-small/src/LibSL/Mesh/VertexFormat_dynamic.cpp
        LibSL-small/src/LibSL/System/System.cpp
        LibSL-small/src/LibSL/Image/Image.cpp
        LibSL-small/src/LibSL/CppHelpers/CppHelpers.cpp
        LibSL-small/src/LibSL/StlHelpers/StlHelpers.cpp
        LibSL-small/src/loki/NullType.h
        LibSL-small/src/loki/static_check.h
        LibSL-small/src/loki/Typelist.h
        LibSL-small/src/loki/TypelistMacros.h
        LibSL-small/src/loki/TypeManip.h
        )

#if (CURVI_BUILD_WITH_OSPQ)
#ADD_EXECUTABLE(curvislice_osqp
#        main.cpp
#        gcode.h
#        gcode.cpp
#        helpers.h
#        TetMesh.h
#        TetMesh.cpp
#        thicknesses.h
#        MeshFormat_msh.h
#        MeshFormat_msh.cpp
#        ${LIBSL_SOURCES}
#        )
#endif()

ADD_EXECUTABLE(curvislice_grb
        main.cpp
        gcode.h
        gcode.cpp
        helpers.h
        TetMesh.h
        TetMesh.cpp
        thicknesses.h
        MeshFormat_msh.h
        MeshFormat_msh.cpp
        ${LIBSL_SOURCES}
        )


ADD_EXECUTABLE(uncurve
        uncurve.cpp
        gcode.h
        gcode.cpp
        TetMesh.h
        TetMesh.cpp
        thicknesses.h
        MeshFormat_msh.h
        MeshFormat_msh.cpp
        ${LIBSL_SOURCES}
        )

ADD_EXECUTABLE(reorient
        reorient.cpp
        ${LIBSL_SOURCES}
        )

ADD_EXECUTABLE(subdiv
        subdiv.cpp
        ${LIBSL_SOURCES}
        )

#if (CURVI_BUILD_WITH_OSPQ)
#target_compile_options(curvislice_osqp PRIVATE -DOSQP)
#endif()

target_compile_options(curvislice_grb PRIVATE -DGRB)
target_compile_options(uncurve PRIVATE -DGRB)

if(WIN32)
#if (CURVI_BUILD_WITH_OSPQ)
#    target_link_libraries(curvislice_osqp shlwapi osqpstatic SolverWrapper)
#    target_link_libraries(curvislice_osqp debug gurobi_c++mdd2017)
#    target_link_libraries(curvislice_osqp optimized gurobi_c++md2017)
#endif()
    target_link_libraries(curvislice_grb shlwapi gurobi90)
    target_link_libraries(curvislice_grb debug gurobi_c++mdd2017)
    target_link_libraries(curvislice_grb optimized gurobi_c++md2017)
    target_link_libraries(reorient shlwapi)
    target_link_libraries(subdiv shlwapi)
    target_link_libraries(uncurve shlwapi gurobi90)
    target_link_libraries(uncurve debug gurobi_c++mdd2017)
    target_link_libraries(uncurve optimized gurobi_c++md2017)
else()
#if (CURVI_BUILD_WITH_OSPQ)
#    target_link_libraries(curvislice_osqp osqpstatic SolverWrapper)
#    target_link_libraries(curvislice_osqp -L/$ENV{GUROBI_HOME}/lib/ -lgurobi_c++ -lgurobi90)
#endif()
    target_link_libraries(curvislice_grb -L/$ENV{GUROBI_HOME}/lib/ -lgurobi_c++ -lgurobi90)
    target_link_libraries(uncurve -L/$ENV{GUROBI_HOME}/lib/ -lgurobi_c++ -lgurobi90)
endif(WIN32)

#if (CURVI_BUILD_WITH_OSPQ)
#set_target_properties(curvislice_osqp PROPERTIES
#        CXX_STANDARD 17
#        CXX_EXTENSIONS OFF
#        )
#endif()

set_target_properties(curvislice_grb PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

set_target_properties(uncurve PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

set_target_properties(reorient PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

#if (CURVI_BUILD_WITH_OSPQ)
#INSTALL(TARGETS curvislice_osqp
#        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
#        )
#endif()

INSTALL(TARGETS curvislice_grb
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )

INSTALL(TARGETS reorient
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )

INSTALL(TARGETS subdiv
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )

INSTALL(TARGETS uncurve
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )
