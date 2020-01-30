CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ADD_EXECUTABLE( TetViz
	TetViz.h
	TetViz.cpp
	shade.fp
	shade.vp
	shade.h
	slicerror.fp
	slicerror.vp
	slicerror.h
	gcode.h
	gcode.cpp
	TetMesh.h
	TetMesh.cpp  
)

TARGET_LINK_LIBRARIES(TetViz
	LibSL
	LibSL_gl4
)

AUTO_BIND_SHADERS(
	shade
  slicerror
)

INSTALL(TARGETS TetViz
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}
)
