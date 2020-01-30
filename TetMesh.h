#pragma once

#include <vector>
#include <algorithm>

#include <LibSL/LibSL.h>
#include <LibSL/Mesh/Mesh.h>
#include <LibSL/Math/Tuple.h>

#include "MeshFormat_msh.h"

using namespace LibSL::Math;
using namespace std;


typedef Tuple<double, 9> M3x3;

class TetMesh
{
private :
  AutoPtr<MeshFormat_msh> m_mesh;

  vector<pair<v4u, bool>>  m_tetrahedrons;
  vector<v3f>  m_vertices;
  vector<v3u>  m_triangles;
  vector<uint> m_surfaces;

  map<uint, vector<uint>> m_tet_neighborhood;
  map<uint, vector<uint>> m_tri_neighborhood;
  map<uint, vector<uint>> m_ver_to_tris;
public:
  TetMesh();
  TetMesh(TetMesh*);
  ~TetMesh();

  static TetMesh* load(const char *fname);
  static void save(const char *fname, TetMesh *mesh);

  size_t numVertices();
  v3f& vertexAt(uint n);

  size_t numTetrahedrons();
  bool isInside(uint t);
  v4u& tetrahedronAt(uint t);
  vector<uint>& tet_neighbours(uint t);
  vector<uint>& tri_neighbours(uint t);
  vector<uint>  ver_neighbours(uint t);

  uint getTetrahedronSurface(uint t);
  
  size_t numTriangles();
  v3u& triangleAt(uint t);

  size_t numSurfaces();
  v3u& surfaceTriangleAt(uint s);


  M3x3  getGradientMatrix(uint t);
  AABox getBBox();

  void reorientSurface();

};
