#include "TetMesh.h"

//---------------------------------------------------------------------------
#include "LibSL.precompiled.h"
//---------------------------------------------------------------------------

using namespace LibSL::Mesh;

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;
#include <LibSL/Math/Vertex.h>
using namespace LibSL::Math;
#include <LibSL/CppHelpers/BasicParser.h>
using namespace LibSL::CppHelpers;

#include <string>
#include <algorithm>

#include <set>
#include <queue>



TetMesh::TetMesh()
{
}

TetMesh::TetMesh(TetMesh* mesh) 
{

  for (auto vertex : mesh->m_vertices) {
    m_vertices.emplace_back(vertex);
  }
  for (auto tet : mesh->m_tetrahedrons) {
    m_tetrahedrons.emplace_back(tet);
  }
  for (auto tri : mesh->m_triangles) {
    m_triangles.emplace_back(tri);
  }
  for (auto srfs : mesh->m_surfaces) {
    m_surfaces.emplace_back(srfs);
  }

}


TetMesh::~TetMesh()
{
}

size_t TetMesh::numVertices() {
  return m_vertices.size();
}

v3f& TetMesh::vertexAt(uint n) {
  return m_vertices.at(n);
}

size_t TetMesh::numTetrahedrons() {
  return m_tetrahedrons.size();
}

bool TetMesh::isInside(uint t) {
  return !m_tetrahedrons.at(t).second;
}

v4u& TetMesh::tetrahedronAt(uint t) {
  return m_tetrahedrons.at(t).first;
}

size_t TetMesh::numTriangles() {
  return m_triangles.size();
}

v3u& TetMesh::triangleAt(uint t) {
  return m_triangles.at(t);
}

size_t TetMesh::numSurfaces() {
  return m_surfaces.size();
}

v3u& TetMesh::surfaceTriangleAt(uint s) {
  return m_triangles.at(m_surfaces.at(s));
}


using namespace std;

//---------------------------------------------------------------------------



TetMesh* TetMesh::load(const char *fname)
{
  TetMesh* mesh = new TetMesh();

  mesh->m_mesh = AutoPtr<MeshFormat_msh>(new MeshFormat_msh);
  mesh->m_mesh->parse(fname);

  MeshFormat_msh::msh_format_info format = mesh->m_mesh->m_format;
  MeshFormat_msh::msh_entity_data element_inout;


  // Check if the domain is tetrahedralized or just the object
  bool domain = false;
  for (auto data : format.element_data) {
    if (data.string_tags[0] == "\"in/out\"") {
      element_inout = data;
      domain = true;
    }
  }
   
  map<v3u, pair<v3u, uint>> tris;
  vector<uint> srfTris;

  int nb_elements = format.elements.size();

  srfTris.resize(nb_elements * 4); // assuming that the file only contains tets

  // faces
  uint fi = 0;
  ForArray(format.elements, n_tet) { // TODO: manage non tet
    bool in = domain && element_inout.values[n_tet][0] == 1.f;

    mesh->m_tetrahedrons.push_back(make_pair(
      v4u(format.elements.at(n_tet)[0] - 1,
          format.elements.at(n_tet)[1] - 1,
          format.elements.at(n_tet)[2] - 1,
          format.elements.at(n_tet)[3] - 1),
      in));

    ForIndex(off, 4) {
      uint A = format.elements.at(n_tet)[(off + 0) % 4] - 1;
      uint B = format.elements.at(n_tet)[(off + 1) % 4] - 1;
      uint C = format.elements.at(n_tet)[(off + 2) % 4] - 1;

      uint _A = A, _B = B, _C = C;

      if (A > B) std::swap(A, B);
      if (B > C) std::swap(B, C);
      if (A > B) std::swap(A, B);

      // extract the surface
      v3u t(A, B, C);
      v3u r(_A, _B, _C);

      auto it = tris.find(t);
      if (it == tris.end()) {
        tris[t].first = r;
        if (in) {
          tris[t].second = 2;
        } else {
          tris[t].second = 1;
        }
        fi++;
      } else {
        if (in) {
          tris[t].first   = r;
          tris[t].second += 2;
        } else {
          tris[t].second += 1;
        }
      }
    }
  }

  //vertices
  ForArray(mesh->m_mesh->m_format.nodes, n_ver) {
    v3f pos = mesh->m_mesh->m_format.nodes.at(n_ver);
    mesh->m_vertices.push_back(pos);
  }

  // tris & srf tris
  uint i = 0;
  for (auto it : tris) {
    mesh->m_triangles.push_back(it.second.first);
    if (it.second.second % 2)
      mesh->m_surfaces.push_back(i);
    i++;
  }

  mesh->reorientSurface();

  // done
  return (mesh);
}



void TetMesh::save(const char *fname, TetMesh *mesh) {
  FILE *f = NULL;
  fopen_s(&f, fname, "wb");
  if (f == NULL) {
    throw Fatal("[MeshFormat_msh::save] - cannot open file '%s' for writing", fname);
  }


  // header
  char header[80];
  memset(header, 0x00, 80);
  sprintf(header, "LibSL_STL_1.0");
  fwrite(header, sizeof(char), 80, f);
  // num triangles
  size_t nTris = mesh->numSurfaces();
  fwrite(&nTris, sizeof(uint), 1, f);
  // vertices
  ForIndex(t, nTris) {
    v3u tri = mesh->surfaceTriangleAt(t);
    v3f nrm = cross(mesh->vertexAt(tri[1]) - mesh->vertexAt(tri[0]), mesh->vertexAt(tri[2]) - mesh->vertexAt(tri[0]));
    nrm = normalize_safe(nrm);
    v3f pt[3];
    ForIndex(v, 3) {
      pt[v] = mesh->vertexAt(tri[v]);
    }
    fwrite(&nrm, sizeof(float), 3, f);
    fwrite(&pt, sizeof(float), 3 * 3, f);
    ushort attr = 0;
    fwrite(&attr, sizeof(ushort), 1, f);
  }

  fclose(f);


}















/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
int LUPDecompose(double **A, int N, double Tol, int *P) {

  int i, j, k, imax;
  double maxA, absA;

  for (i = 0; i <= N; i++)
    P[i] = i; //Unit permutation matrix, P[N] initialized with N

  for (i = 0; i < N; i++) {
    maxA = 0.0;
    imax = i;

    for (k = i; k < N; k++)
      if ((absA = fabs(A[k][i])) > maxA) {
        maxA = absA;
        imax = k;
      }

    if (maxA < Tol) return 0; //failure, matrix is degenerate

    if (imax != i) {
      //pivoting P
      swap(P[i], P[imax]);

      //pivoting rows of A
      swap(A[i], A[imax]);

      //counting pivots starting from N (for determinant)
      P[N]++;
    }

    for (j = i + 1; j < N; j++) {
      A[j][i] /= A[i][i];

      for (k = i + 1; k < N; k++)
        A[j][k] -= A[j][i] * A[i][k];
    }
  }

  return 1;  //decomposition done 
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
* OUTPUT: x - solution vector of A*x=b
*/
void LUPSolve(double **A, int *P, double *b, int N, double *x) {

  for (int i = 0; i < N; i++) {
    x[i] = b[P[i]];

    for (int k = 0; k < i; k++)
      x[i] -= A[i][k] * x[k];
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++)
      x[i] -= A[i][k] * x[k];

    x[i] = x[i] / A[i][i];
  }
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
* OUTPUT: IA is the inverse of the initial matrix
*/
void LUPInvert(double **A, int *P, int N, double **IA) {

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      if (P[i] == j)
        IA[i][j] = 1.0;
      else
        IA[i][j] = 0.0;

      for (int k = 0; k < i; k++)
        IA[i][j] -= A[i][k] * IA[k][j];
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int k = i + 1; k < N; k++)
        IA[i][j] -= A[i][k] * IA[k][j];

      IA[i][j] = IA[i][j] / A[i][i];
    }
  }
}

/* INPUT: A,P filled in LUPDecompose; N - dimension.
* OUTPUT: Function returns the determinant of the initial matrix
*/
double LUPDeterminant(double **A, int *P, int N) {

  double det = A[0][0];

  for (int i = 1; i < N; i++)
    det *= A[i][i];

  if ((P[N] - N) % 2 == 0)
    return det;
  else
    return -det;
}

#include "gurobi_c++.h"
M3x3 TetMesh::getGradientMatrix(uint t) {
  M3x3 m(0.);

  v3f A = vertexAt(tetrahedronAt(t)[0]);
  v3f B = vertexAt(tetrahedronAt(t)[1]);
  v3f C = vertexAt(tetrahedronAt(t)[2]);
  v3f D = vertexAt(tetrahedronAt(t)[3]);

  A -= D;
  B -= D;
  C -= D;

  static GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_OutputFlag, false);

  GRBVar M[9];

  ForIndex(i, 9) {
    M[i] = model.addVar(-1000.0f, 1000.0f, m[i], GRB_CONTINUOUS);
  }

  model.addConstr(M[0] * A[0] + M[1] * B[0] + M[2] * C[0] == 1.);
  model.addConstr(M[0] * A[1] + M[1] * B[1] + M[2] * C[1] == 0.);
  model.addConstr(M[0] * A[2] + M[1] * B[2] + M[2] * C[2] == 0.);

  model.addConstr(M[3] * A[0] + M[4] * B[0] + M[5] * C[0] == 0.);
  model.addConstr(M[3] * A[1] + M[4] * B[1] + M[5] * C[1] == 1.);
  model.addConstr(M[3] * A[2] + M[4] * B[2] + M[5] * C[2] == 0.);

  model.addConstr(M[6] * A[0] + M[7] * B[0] + M[8] * C[0] == 0.);
  model.addConstr(M[6] * A[1] + M[7] * B[1] + M[8] * C[1] == 0.);
  model.addConstr(M[6] * A[2] + M[7] * B[2] + M[8] * C[2] == 1.);

  model.optimize();
  
  ForIndex(i, 9) {
    m[i] = M[i].get(GRB_DoubleAttr_X);
  }
  


  // without optim
  
  /*
  double** M_2 = new double*[9];
  for (int i = 0; i < 3; ++i) {
    M_2[i] = new double[9];
  }

  double b[9] = { 1.0 , 0.0 , 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  double res[9] = { 0, 0, 0 , 0, 0, 0, 0, 0, 0};

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) M_2[i][0] = A[i];
    for (int i = 0; i < 3; i++) M_2[i][1] = B[i];
    for (int i = 0; i < 3; i++) M_2[i][2] = C[i];

    int permut[4] = { 0, 0, 0, 0 };
    LUPDecompose(M_2, 3, 1E-18, permut);

    
    LUPSolve(M_2, permut, &b[k * 3], 3, &res[k * 3]);
  }

  for (int i = 0; i < 3; ++i) {
    delete[] M_2[i];
  }
  delete[] M_2;

  for(int i = 0; i < 9; i++)
    m[i] = res[i];
  */

  return m;
}

vector<uint>& TetMesh::tet_neighbours(uint t) {
  if (m_tet_neighborhood.empty()) {
    for(int i = 0; i < m_tetrahedrons.size()-1; i++) {
      v4u origin = tetrahedronAt(i);
      for (int j = i + 1; j < m_tetrahedrons.size(); j++) {
        v4u curr = tetrahedronAt(j);
        int common = 0;
        ForIndex(a, 4) {
          ForIndex(b, 4) {
            if (curr[a] == origin[b]) {
              common++;
            }
          }
        }
        if (common == 3) {
          m_tet_neighborhood[i].push_back(j);
          m_tet_neighborhood[j].push_back(i);
        }
      }
    }
  }

  return m_tet_neighborhood[t];
}

vector<uint>& TetMesh::tri_neighbours(uint t) {
  if (m_tri_neighborhood.empty()) {
    for (int i = 0; i < numSurfaces() - 1; i++) {
      v3u origin = surfaceTriangleAt(i);
      for (int j = i + 1; j < numSurfaces(); j++) {
        v3u curr = surfaceTriangleAt(j);
        int common = 0;
        ForIndex(a, 3) {
          ForIndex(b, 3) {
            if (curr[a] == origin[b]) {
              common++;
            }
          }
        }
        if (common == 2) {
          m_tri_neighborhood[i].push_back(j);
          m_tri_neighborhood[j].push_back(i);
        }
      }
    }
  }

  return m_tri_neighborhood[t];
}

vector<uint> TetMesh::ver_neighbours(uint t) {
  if (m_ver_to_tris.empty()) {
    ForIndex (s, m_surfaces.size()) {
      ForIndex(v, 3) {
        m_ver_to_tris[v].push_back(s);
      }
    }
  }

  vector<uint> vec;

  for (uint s : m_ver_to_tris[t]) {
    ForIndex(v, 3) {
      if (m_triangles.at(s)[v] != t) {
        vec.push_back(m_triangles.at(s)[v]);
      }
    }
  }

  sort(vec.begin(), vec.end());
  vec.erase(unique(vec.begin(), vec.end()), vec.end());

  return vec;
}

uint TetMesh::getTetrahedronSurface(uint t) {
  v3u tri = m_triangles.at(t);
  for (int i = 0; i < m_tetrahedrons.size(); i++) {
    if (!m_tetrahedrons.at(i).second) continue;
    v4u tet = m_tetrahedrons.at(i).first;
    int common = 0;
    ForIndex(a, 4) {
      ForIndex(b, 3) {
        if (tet[a] == tri[b]) {
          common++;
        }
      }
    }
    if (common == 3) {
      return i;
    }
  }
  return 0;
}

AABox TetMesh::getBBox() {
  AABox bbox;
  for (v3f vertex : m_vertices) {
    bbox.addPoint(vertex);
  }
  return bbox;
}

void TetMesh::reorientSurface()
{
  if (numSurfaces() == 0) return;
  Array<bool> visited(numSurfaces());
  visited.fill(false);

  // build edge info
  map<v2i, vector<int> > edge_to_tris;
  ForIndex(t, numSurfaces()) {
    v3u tri = surfaceTriangleAt(t);
    ForIndex(i, 3) {
      v2i e = v2i(min(tri[i], tri[(i + 1) % 3]), max(tri[i], tri[(i + 1) % 3]));
      edge_to_tris[e].push_back(t);
    }
  }

  while (true) {

    // find a lowest triangle not visited
    int next = -1;
    float minz = FLT_MAX;
    ForIndex(t, numSurfaces()) {
      if (!visited[t]) {
        ForIndex(i, 3) {
          float z = vertexAt(surfaceTriangleAt(t)[i])[2];
          if (z < minz) {
            minz = z;
            next = t;
          }
        }
      }
    }
    if (next == -1) {
      break; // done!
    }

    // orient 'next'
    v3u tn = surfaceTriangleAt(next);
    v3f nrm = cross(vertexAt(tn[1]) - vertexAt(tn[0]), vertexAt(tn[2]) - vertexAt(tn[0]));
    if (nrm[2] > 0) {
      swap(tn[0], tn[2]);
      surfaceTriangleAt(next) = tn;
    }

    // orient by local growth
    std::queue<int> q;
    q.push(next);
    visited[next] = true;

    while (!q.empty()) {
      // pop current
      int curid = q.front();
      v3u cur = surfaceTriangleAt(curid);
      q.pop();
      // get neighbors
      ForIndex(i, 3) {
        v2i e = v2i(min(cur[i], cur[(i + 1) % 3]), max(cur[i], cur[(i + 1) % 3]));
        const vector<int>& neighs = edge_to_tris[e];
        if (neighs.size() == 2) { // only trust two-manifold edges
          ForIndex(n, neighs.size()) {
            if (!visited[neighs[n]]) {
              visited[neighs[n]] = true;
              // check orientation
              v3u tri_n = surfaceTriangleAt(neighs[n]);
              bool fliped = false;
              v2i e_cur = v2i(cur[i], cur[(i + 1) % 3]);
              ForIndex(j, 3) {
                v2i e_j = v2i(tri_n[j], tri_n[(j + 1) % 3]);
                if (e_j == e_cur) {
                  fliped = true; break;
                }
              }
              if (fliped) {
                // flip
                swap(tri_n[0], tri_n[2]);
                surfaceTriangleAt(neighs[n]) = tri_n;
              }
              q.push(neighs[n]);
            }
          }
        }
      }
    }
  }

  // search for face at bottom to check if we have to flip all
  float bottomz = std::numeric_limits<float>::max();
  ForIndex(t, numSurfaces()) {
    ForIndex(i, 3) {
      v3f pt = vertexAt(surfaceTriangleAt(t)[i]);
      bottomz = min(bottomz, pt[2]);
    }
  }
  bool flip = false;
  ForIndex(t, numSurfaces()) {
    bool bottom = true;
    v3f pts[3];
    ForIndex(i, 3) {
      v3f pt = vertexAt(surfaceTriangleAt(t)[i]);
      if (abs(pt[2] - bottomz) > 1e-3f) {
        bottom = false;
      }
      pts[i] = pt;
    }
    if (bottom) {
      // compute normal to check it faces down
      v3f nrm = cross(pts[1] - pts[0], pts[2] - pts[0]);
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
    ForIndex(t, numSurfaces()) {
      v3u tri = surfaceTriangleAt(t);
      std::swap(tri[0], tri[1]);
      surfaceTriangleAt(t) = tri;
    }
  }
}