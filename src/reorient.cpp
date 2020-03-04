#include <algorithm>
#include <LibSL/LibSL.h>
#include <LibSL/Memory/ArrayTools.h>
#include <LibSL/Mesh/MeshFormat_stl.h>

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  sl_assert(argc == 2);
  TriangleMesh_Ptr mesh( loadTriangleMesh(argv[1]) );
  // merge all vertices
  mesh->mergeVerticesExact();
  // reorient
  mesh->reorientTriangles();
  // search for face at bottom to check if we have to flip all
  AAB<3> bx = mesh->bbox();
  bool flip = false;
  ForIndex(t, mesh->numTriangles()) {
    bool bottom = true;
    v3d pts[3];
    ForIndex(i, 3) { 
      v3d pt = v3d(mesh->posAt(mesh->triangleAt(t)[i])); 
      if (abs(pt[2] - bx.minCorner()[2]) > 1e-3f) {
        bottom = false;
      }    
      pts[i] = pt;
    }
    if (bottom) {
      // compute normal to check it faces down
      v3d nrm = cross(pts[1] - pts[0], pts[2] - pts[0]);
      if (length(nrm) > 1e-6f) {
        // can trust normal
        if (nrm[2] > 0) {
          // flip all
          flip = true;
        }
        // done
        break;
      }
    }
  }
  if (flip) {
    ForIndex(t, mesh->numTriangles()) {
      v3u tri = mesh->triangleAt(t);
      std::swap(tri[0], tri[1]);
      mesh->triangleAt(t) = tri;
    }
  }
  // save
  // Not declared on windows (LIBSL/System/System.cpp)
  //CreateDirectoryA("scratch", NULL);
  saveTriangleMesh("scratch\\oriented.stl",mesh.raw());
  return 0;
}
