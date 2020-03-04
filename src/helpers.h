
#include <LibSL/LibSL.h>
// --------------------------------------------------------------
// --------------    FUNCTIONS / HELPERS   ----------------------
// --------------------------------------------------------------

// ScTP computes the scalar triple product
#define ScTP(a, b, c) dot(a, cross(b,c))

// --------------------------------------------------------------

float tetrahedron_signed_volume(const Tuple<v3f, 4>& tet)
{
  v3f vab = tet[1] - tet[0];
  v3f vac = tet[2] - tet[0];
  v3f vad = tet[3] - tet[0];

  return 1.0F / 6.0F * ScTP(vab, vac, vad);
}

// --------------------------------------------------------------

float tetrahedron_signed_volume(TetMesh& mesh, const v4u& tet)
{
  return tetrahedron_signed_volume(Tuple<v3f, 4>(mesh.vertexAt(tet[0]), mesh.vertexAt(tet[1]), mesh.vertexAt(tet[2]), mesh.vertexAt(tet[3])));
}

// --------------------------------------------------------------

float tetrahedron_volume(const Tuple<v3f, 4>& tet)
{
  return abs(tetrahedron_signed_volume(tet));
}

// --------------------------------------------------------------

float tetrahedron_volume(TetMesh& mesh, const v4u& tet)
{
  return tetrahedron_volume(Tuple<v3f, 4>(mesh.vertexAt(tet[0]), mesh.vertexAt(tet[1]), mesh.vertexAt(tet[2]), mesh.vertexAt(tet[3])));
}

// --------------------------------------------------------------

v4f tetrahedron_barycenter_coefs(const Tuple<v3f, 4>& tet, v3f pt)
{
  v3f vap = pt - tet[0];
  v3f vbp = pt - tet[1];

  v3f vab = tet[1] - tet[0];
  v3f vac = tet[2] - tet[0];
  v3f vad = tet[3] - tet[0];

  v3f vbc = tet[2] - tet[1];
  v3f vbd = tet[3] - tet[1];

  float va = 1.0F / 6.0F * ScTP(vbp, vbd, vbc);
  float vb = 1.0F / 6.0F * ScTP(vap, vac, vad);
  float vc = 1.0F / 6.0F * ScTP(vap, vad, vab);
  float vd = 1.0F / 6.0F * ScTP(vap, vab, vac);

  return v4f(va, vb, vc, vd) / tetrahedron_volume(tet);
}

// --------------------------------------------------------------

float triangle_area(const Tuple<v3d, 3>& tri)
{
  v3d vab = tri[1] - tri[0];
  v3d vac = tri[2] - tri[0];

  float v = 1.0F / 2.0F * length(cross(vab, vac));

  return abs(v);
}

// --------------------------------------------------------------

float triangle_area(TetMesh& mesh, const v3u& tri)
{
  return triangle_area(Tuple<v3d, 3>(v3d(mesh.vertexAt(tri[0])), v3d(mesh.vertexAt(tri[1])), v3d(mesh.vertexAt(tri[2]))));
}

// --------------------------------------------------------------

float triangle_area_xy(const Tuple<v3d, 3>& tri)
{
  v3d vab = tri[1] - tri[0];
  v3d vac = tri[2] - tri[0];
  vab[2] = 0.0F; vac[2] = 0.0F;
  float v = 1.0F / 2.0F * length(cross(vab, vac));
  return abs(v);
}

// --------------------------------------------------------------

double triangle_area_xy(TetMesh& mesh, const v3u& tri)
{
  return triangle_area_xy(Tuple<v3d, 3>(v3d(mesh.vertexAt(tri[0])), v3d(mesh.vertexAt(tri[1])), v3d(mesh.vertexAt(tri[2]))));
}

// --------------------------------------------------------------

float component_area(TetMesh& mesh, const vector<int>& surfaces)
{
  float area = 0.0F;
  for (int surface : surfaces) {
    area += triangle_area(mesh, mesh.surfaceTriangleAt(surface));
  }
  return area;
}

// --------------------------------------------------------------

template<typename T_var>
Tuple<SLRExpr<double>, 3> tetrahedron_gradient(TetMesh& mesh, map<uint, T_var>& h, const Array<M3x3>& mats, uint t)
{
  v4u tet = mesh.tetrahedronAt(t);

  T_var ha = h[tet[0]];
  T_var hb = h[tet[1]];
  T_var hc = h[tet[2]];
  T_var hd = h[tet[3]];

    SLRExpr<double> dhdx = (ha - hd) * mats[t][0]
                  + (hb - hd) * mats[t][1]
                  + (hc - hd) * mats[t][2];

    SLRExpr<double> dhdy = (ha - hd) * mats[t][3]
                  + (hb - hd) * mats[t][4]
                  + (hc - hd) * mats[t][5];

    SLRExpr<double> dhdz = (ha - hd) * mats[t][6]
                  + (hb - hd) * mats[t][7]
                  + (hc - hd) * mats[t][8];

  return Tuple<SLRExpr<double>, 3>(dhdx, dhdy, dhdz);
}

// --------------------------------------------------------------

template < class T_Graph >
class TestFlattened
{
private:
  const set<int>& flats;
public:
  TestFlattened(const set<int>& flats_) : flats(flats_) {}
  bool operator()(const T_Graph& graph, LibSL::DataStructures::t_NodeId node) const
  {
    return (flats.find(node) != flats.end());
  }
};

// ---------------------------------------------------------------

void connected_components(TetMesh& mesh, const set<int>& surfaces_to_flatten, vector<vector<int> >& _comps)
{
  typedef LibSL::DataStructures::Graph<int, Loki::NullType, false> t_Graph;
  t_Graph g;
  ForIndex(v, mesh.numSurfaces()) {
    g.addNode(v);
  }
  ForIndex(v, mesh.numSurfaces()) {
    vector<uint> neighs = mesh.tri_neighbours(v);
    for (auto n : neighs) {
      g.addEdge(v, n, Loki::NullType());
    }
  }
  LibSL::DataStructures::GraphAlgorithms::ConnectedComponents<t_Graph, LibSL::DataStructures::GraphAlgorithms::DefaultEdgeTester<t_Graph>, TestFlattened<t_Graph>> comps;
  Array<bool> visited;
  vector<LibSL::DataStructures::t_NodeId> connected;
  _comps.clear();
  while (comps.enumerateConnectedComponents(g, visited, connected,
         LibSL::DataStructures::GraphAlgorithms::DefaultEdgeTester<t_Graph>(),
         TestFlattened<t_Graph>(surfaces_to_flatten))) {
    _comps.push_back(vector<int>());
    for (auto n : connected) {
      _comps.back().push_back(n);
    }
  }
}

// ---------------------------------------------------------------

double triangle_non_flatness(v3u tri, SLRModel<double> *model)
{
  double va = model->getVarByName(sprint("z_%03d", tri[0])).get();
  double vb = model->getVarByName(sprint("z_%03d", tri[1])).get();
  double vc = model->getVarByName(sprint("z_%03d", tri[2])).get();
  double non_flatness = max((va - vb)*(va - vb), max((va - vc)*(va - vc), (vb - vc)*(vb - vc)));
  return non_flatness;
}

// --------------------------------------------------------------

double component_non_flatness(TetMesh& mesh, const std::vector<int>& tris, SLRModel<double> *model, double max_admissible)
{
  double non_flatness = 0.0, max_non_flatness = 0.0;
  double comp_area = 0.0;
  for (auto t : tris) {
    v3u tri = mesh.surfaceTriangleAt(t);
    double ta = triangle_area_xy(mesh, tri);
    double nftri = triangle_non_flatness(tri, model);
    max_non_flatness = max(max_non_flatness, nftri);
    non_flatness += nftri * ta;
    comp_area += ta;
  }
  if (max_non_flatness > max_admissible) {
    return max_non_flatness;
  } else {
    return non_flatness / comp_area;
  }
}
