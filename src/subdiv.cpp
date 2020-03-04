#include <algorithm>
#include <LibSL/LibSL.h>
#include <LibSL/Memory/ArrayTools.h>
#include <LibSL/Mesh/MeshFormat_stl.h>

#include <iostream>

using namespace std;

// --------------------------------------------------------------

void subdivide(
  v3f p0, v3f p1, v3f p2,
  float minlen,
  vector<v3f>& _pts, vector<v3u>& _tris)
{
  v3f pts[3] = { p0,p1,p2 };
  float sqm = (minlen) * (minlen);
  ForIndex(i, 3) {
    int a = i;
    int b = (i + 1) % 3;
    int c = (i + 2) % 3;
    if (sqLength(pts[a] - pts[b]) > sqm) {
      // split
      //         a
      //         /\
      //        /  \
      //   abh /    \
      //      / \    \
      //     /     \  \
      //  b /---------\\ c
      //        
      v3f abh = pts[a] + (pts[b] - pts[a]) / 2.0f;
      subdivide(pts[c], pts[a], abh, minlen, _pts, _tris);
      subdivide(pts[b], pts[c], abh, minlen, _pts, _tris);
      return;
    }
  }
  int id = (int)_pts.size();
  _pts.push_back(p0);
  _pts.push_back(p1);
  _pts.push_back(p2);
  _tris.push_back(v3u(id, id + 1, id + 2));
}

// --------------------------------------------------------------

int main(int argc, char **argv)
{
  sl_assert(argc == 2);
  TriangleMesh_Ptr mesh( loadTriangleMesh(argv[1]) );
  // merge all vertices
  mesh->mergeVerticesExact();
  // reorient
  vector<v3f> vpts;
  vector<v3u> vtris;
  ForIndex(t, mesh->numTriangles()) {
    v3u tri = mesh->triangleAt(t);
    v3f p0 = mesh->posAt(tri[0]);
    v3f p1 = mesh->posAt(tri[1]);
    v3f p2 = mesh->posAt(tri[2]);
    subdivide(p0, p1, p2, 1.0f/*mm*/, vpts, vtris);
  }
  TriangleMesh_Ptr tmesh(new TriangleMesh_generic<MeshFormat_stl::t_VertexData>(vpts.size(), vtris.size()));
  ForIndex(v, tmesh->numVertices()) {
    tmesh->posAt(v) = v3f(vpts[v]);
  }
  ForIndex(t, tmesh->numTriangles()) {
    tmesh->triangleAt(t) = vtris[t];
  }
  // save
  // Not declared on windows (LIBSL/System/System.cpp)
  //CreateDirectoryA("scratch",NULL);
  saveTriangleMesh("scratch\\subdiv.stl", tmesh.raw());
  return 0;
}

// --------------------------------------------------------------
