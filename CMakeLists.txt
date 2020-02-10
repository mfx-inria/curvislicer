cmake_minimum_required(VERSION 2.6)
project(curvislice)

OPTION(CURVI_BUILD_WITH_GRB  "Build with GRB (stable solver)" OFF)

add_subdirectory(libs/SolverWrapper)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}/libs/LibSL-small/src
        ${PROJECT_SOURCE_DIR}/libs/LibSL-small/src/LibSL
        ${PROJECT_SOURCE_DIR}/libs/tclap/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/include/
        ${PROJECT_SOURCE_DIR}/libs/SolverWrapper/src/
)

CONFIGURE_FILE(
        "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
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

ADD_EXECUTABLE(curvislice_osqp
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

if (CURVI_BUILD_WITH_GRB)
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
endif(CURVI_BUILD_WITH_GRB)


ADD_EXECUTABLE(reorient
        reorient.cpp
        ${LIBSL_SOURCES}
        )

ADD_EXECUTABLE(subdiv
        subdiv.cpp
        ${LIBSL_SOURCES}
        )

target_compile_options(curvislice_osqp PRIVATE -DOSQP)

if (CURVI_BUILD_WITH_GRB)
    target_compile_options(curvislice_grb PRIVATE -DGRB)
    target_compile_options(uncurve PRIVATE -DGRB)
endif(CURVI_BUILD_WITH_GRB)


if(WIN32)
    target_link_libraries(curvislice_osqp shlwapi osqpstatic SolverWrapper)

    target_link_libraries(curvislice_grb shlwapi gurobi81)
    target_link_libraries(curvislice_grb debug gurobi_c++mdd2017)
    target_link_libraries(curvislice_grb optimized gurobi_c++md2017)

    target_link_libraries(reorient shlwapi)
    target_link_libraries(subdiv shlwapi)
    target_link_libraries(uncurve shlwapi gurobi81)
    target_link_libraries(uncurve debug gurobi_c++mdd2017)
    target_link_libraries(uncurve optimized gurobi_c++md2017)
else()
    target_link_libraries(curvislice_osqp SolverWrapper)

    if (CURVI_BUILD_WITH_GRB)
        target_link_libraries(curvislice_grb SolverWrapper)
        target_link_libraries(uncurve SolverWrapper)
    endif(CURVI_BUILD_WITH_GRB)
endif(WIN32)

set_target_properties(curvislice_osqp PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )
if (CURVI_BUILD_WITH_GRB)
    set_target_properties(curvislice_grb PROPERTIES
            CXX_STANDARD 17
            CXX_EXTENSIONS OFF
            )

    set_target_properties(uncurve PROPERTIES
            CXX_STANDARD 17
            CXX_EXTENSIONS OFF
            )
endif(CURVI_BUILD_WITH_GRB)


set_target_properties(reorient PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

INSTALL(TARGETS curvislice_osqp
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )

if (CURVI_BUILD_WITH_GRB)
    INSTALL(TARGETS curvislice_grb
            RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
            )

    INSTALL(TARGETS uncurve
            RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
            )
endif(CURVI_BUILD_WITH_GRB)


INSTALL(TARGETS reorient
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )

INSTALL(TARGETS subdiv
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )


