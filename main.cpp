// --------------------------------------------------------------
// ------------------    INCLUDES    ----------------------------
// --------------------------------------------------------------

#include <algorithm>

#include <LibSL/LibSL.h>
#include <LibSL/Memory/ArrayTools.h>

#include <iostream>
#include <queue>
#include <iomanip>
#include <set>
#include <iterator>

#include "gurobi_c++.h"

#include "path.h"
#include "gcode.h"
#include "TetMesh.h"
#include "MeshFormat_msh.h"

#include "thicknesses.h"
#include "helpers.h"


// --------------------------------------------------------------
// ------------------     GUROBI     ----------------------------
// --------------------------------------------------------------

#define PRESOLVE -1 // -1 = auto, 0 = off, 1 = conservative, 2 = agressive

// --------------------------------------------------------------
// ------------------   PARAMETERS   ----------------------------
// --------------------------------------------------------------

#define WEIGHT_OBJ_SLOPE         0.1
#define WEIGHT_OBJ_SMOOTHNESS    0.02
#define WEIGHT_OBJ_FLATTENING    30.0

#define CSTRT_COLLISION_SLOPE    1
#define CSTRT_THICKNESS          1
#define CSTRT_FOLDOVER           1

#define CSTRT_ANCHOR_POINT       0
#define CSTRT_ANCHOR_LAYER       1 & !CSTRT_ANCHOR_POINT

#define OBJ_SURFACE_SLOPE        1

#define OBJ_SMOOTHNESS_L1        0
#define OBJ_SMOOTHNESS_L2        1

#define IGNORE_EMPTINESS         0

#define THRESHOLD_AREA_mm2       10.0
#define THRESHOLD_MIN_VOLUME_mm3 0.01

// --------------------------------------------------------------

using namespace std;

// --------------------------------------------------------------
// ------------   COMMANDLINE PARAMETERS   ----------------------
// --------------------------------------------------------------

string folder;
string filename;

double max_theta; // theta max : maximum angle printable
double obj_angle; // phi : objectif angle if not flat or vertical ---- 0: vertical
float normal_threshold = 0.95f;

int target_num_layers = -1;
float layer_thickness = default_layer_thickness;

float compute_time = 0.0f; // time for a compute iteration (0 = no limit)
int     nb_threads = 8;

bool  fabricable_emptiness = false;

// --------------------------------------------------------------

void display_conf() {
  cout << endl << Console::yellow << "Normal threshold: " << normal_threshold << endl << endl;
  if (fabricable_emptiness) {
    cout << endl << Console::yellow << "<<Fabricable emptyness (for supports)>>" << endl << endl;
  }

#if CSTRT_FOLDOVER
  cout << Console::green << "CSTRT_FOLDOVER " << endl;
#else
  cout << Console::red << "NO CSTRT_FOLDOVER" << endl;
#endif

#if CSTRT_THICKNESS
  cout << Console::green << "CSTRT_THICKNESS " << min_thickness << " < " << layer_thickness << " < " << max_thickness << endl;
#else
  cout << Console::red << "NO CSTRT_THICKNESS" << endl;
#endif

#if CSTRT_COLLISION_SLOPE
  cout << Console::green << "CSTRT_COLLISION_SLOPE " << -max_theta << " < ? < " << max_theta << endl;
#else
  cout << Console::red << "NO CSTRT_COLLISION_SLOPE" << endl;
#endif

#if CSTRT_ANCHOR_POINT
  cout << Console::yellow << "CSTRT_ANCHOR_POINT ( /!\\ to much freedom)" << endl;
#else
#if CSTRT_ANCHOR_LAYER
  cout << Console::green << "CSTRT_ANCHOR_LAYER" << endl;
#else
  cout << Console::red << "NO CSTRT_ANCHOR" << endl;
#endif
#endif

#if OBJ_SMOOTHNESS_L1
  cout << Console::cyan << "OBJ_SMOOTHNESS_L1" << endl;
#endif

#if OBJ_SMOOTHNESS_L2
  cout << Console::cyan << "OBJ_SMOOTHNESS_L2" << endl;
#endif

#if OBJ_SURFACE_SLOPE
  cout << Console::cyan << "OBJ_SURFACE_SLOPE angle = " << obj_angle << endl;
#endif
  cout << Console::cyan << "TARGET_NUM_LAYERS " << target_num_layers << endl;
  cout << Console::gray << endl;
}


/****************************************************************/
/*********************     OBJECTIVES     ***********************/
/****************************************************************/

////////////////////////////////////////////////////////////////////
//                                                                //
//                          SMOOTHNESS                            //
//                                                                //
////////////////////////////////////////////////////////////////////

template<typename T_obj, typename T_var>
void obj_smoothness_l2(TetMesh& mesh, map<uint, T_var>& h, const Array<M3x3>& mats, T_obj& obj) {
  double volume = 0.0;

  ForIndex(t, mesh.numTetrahedrons()) {
    volume += tetrahedron_volume(mesh, mesh.tetrahedronAt(t));
  }

  ForIndex(t, mesh.numTetrahedrons()) {
#if IGNORE_EMPTINESS
    if (!mesh.isInside(t)) {
      continue;
    }
#endif

    v4u tet = mesh.tetrahedronAt(t);
    auto gradient = tetrahedron_gradient(mesh, h, mats, t);
    auto dhdx = gradient[0];
    auto dhdy = gradient[1];
    auto dhdz = gradient[2];

    for (uint i_tet : mesh.tet_neighbours(t)) {
      if (i_tet == t) continue;

      v4u tet_n = mesh.tetrahedronAt(i_tet);
      auto gradient_n = tetrahedron_gradient(mesh, h, mats, i_tet);
      auto dhdx2 = gradient_n[0];
      auto dhdy2 = gradient_n[1];
      auto dhdz2 = gradient_n[2];

      double vtets = tetrahedron_volume(mesh, tet) + tetrahedron_volume(mesh, tet_n);

      // weight
      double w = WEIGHT_OBJ_SMOOTHNESS * max(vtets, THRESHOLD_MIN_VOLUME_mm3) / volume;
      obj += w * (dhdz - dhdz2) * (dhdz - dhdz2);
      obj += w * (dhdx - dhdx2) * (dhdx - dhdx2);
      obj += w * (dhdy - dhdy2) * (dhdy - dhdy2);
    }
  }
}

////////////////////////////////////////////////////////////////////
//                                                                //
//                              SLOPE                             //
//                                                                //
////////////////////////////////////////////////////////////////////

template<typename T_obj, typename T_var>
void obj_slope(TetMesh& mesh, double tot_surface,map<uint, T_var>& h, T_obj& obj, set<int>& to_flatten)
{
  // determine normals (only z-component)
  Array<v3d> normals;
  normals.allocate(mesh.numSurfaces());
  ForIndex(t, mesh.numSurfaces()) {
    v3u tri = mesh.surfaceTriangleAt(t);
    v3f p[3] = { mesh.vertexAt(tri[0]), mesh.vertexAt(tri[1]), mesh.vertexAt(tri[2]) };
    normals[t] = v3d(normalize_safe(cross(p[2] - p[0], p[1] - p[0])));
  }

  ForIndex(s, mesh.numSurfaces()) {
    // Get the surface normal
    v3d n = normals[s];
    // Ignore surfaces tagged as "to flatten".
    if (to_flatten.find(s) != to_flatten.end()) continue;
    // Ignore vertical ones
    if (abs(n[2]) < 1e-3) continue;

    v3u tri = mesh.surfaceTriangleAt(s);
    v3d p[3] = { v3d(mesh.vertexAt(tri[0])), v3d(mesh.vertexAt(tri[1])), v3d(mesh.vertexAt(tri[2])) };

    v3d t = cross(n, v3d(0, 0, 1));

    t = normalize_safe(t);
    v3d q = normalize_safe(cross(v3d(0, 0, 1), t));
    v3d u = normalize_safe(cross(n, t)); // triangle project as a line along u

    double alpha = obj_angle * M_PI / 180.;
    v3d ut = cos(alpha) * v3d(0, 0, 1) * sign(u[2]) + sin(alpha) * q;

    double Ax = p[0][0];
    double Ay = p[0][1];
    T_var  Az = h[tri[0]];

    double Bx = p[1][0];
    double By = p[1][1];
    T_var  Bz = h[tri[1]];

    double Cx = p[2][0];
    double Cy = p[2][1];
    T_var  Cz = h[tri[2]];

    auto exprab = -Bz + Az + ((Bx - Ax)*u[0] + (By - Ay)*u[1] + (Bz - Az)*u[2]) * ut[2];
    auto exprac = -Cz + Az + ((Cx - Ax)*u[0] + (Cy - Ay)*u[1] + (Cz - Az)*u[2]) * ut[2];
    auto exprbc = -Cz + Bz + ((Cx - Bx)*u[0] + (Cy - By)*u[1] + (Cz - Bz)*u[2]) * ut[2];

    // weight
    double w = WEIGHT_OBJ_SLOPE * triangle_area(mesh, tri) / tot_surface;
    obj += w * exprab * exprab;
    obj += w * exprac * exprac;
    obj += w * exprbc * exprbc;
  }
}

////////////////////////////////////////////////////////////////////
//                                                                //
//                       SOFT FLATTENING                          //
//                                                                //
////////////////////////////////////////////////////////////////////

template<typename T_obj, typename T_var>
void obj_soft_flattening(TetMesh& mesh, double tot_surface, map<uint, T_var>& h, T_obj& obj, set<int>& surfaces_to_flatten) {
  for (auto s : surfaces_to_flatten) {
    v3u tri = mesh.surfaceTriangleAt(s);
    uint a = tri[0], b = tri[1], c = tri[2];
    T_var ha = h[a];
    T_var hb = h[b];
    T_var hc = h[c];

    // weight
    double w = WEIGHT_OBJ_FLATTENING * triangle_area(mesh, tri) / tot_surface;
    obj += w * (ha - hb)*(ha - hb);
    obj += w * (ha - hc)*(ha - hc);
    obj += w * (hb - hc)*(hb - hc);
  }
}



// --------------------------------------------------------------

bool gurobi_opt(
        TetMesh& mesh,
        const Array<Tuple<double, 9>>& mats,
        Array<float>& n_pts
)
{
  double thres_flatteness_mm              = min_thickness / 8.0;
  double max_admissible_non_flatteness_mm = max_thickness / 2.0;

  display_conf();

  Elapsed time{};
  time.elapsed();

  // compute total surface area
  double tot_surface = 0.;
  ForIndex(s, mesh.numSurfaces()) {
    tot_surface += triangle_area(mesh, mesh.surfaceTriangleAt(s));
  }

  // determine normals (only z-component)
  Array<v3d> normals;
  normals.allocate(mesh.numSurfaces());
  float bottomz = std::numeric_limits<float>::max();
  // i_bottom_shape is the lowest variable at the bottom of the shape
  // this is needed for slice alignment
  int   i_bottom_shape = -1;
  ForIndex(t, mesh.numSurfaces()) {
    v3u tri = mesh.surfaceTriangleAt(t);
    v3f p[3] = { mesh.vertexAt(tri[0]), mesh.vertexAt(tri[1]), mesh.vertexAt(tri[2]) };
    normals[t] = v3d(normalize_safe(cross(p[2] - p[0], p[1] - p[0])));
    ForIndex(i, 3) {
      if (p[i][2] < bottomz) {
        bottomz = p[i][2];
        i_bottom_shape = tri[i];
      }
    }
  }

  // find out lowest domain position and max possible height on solution
  int min_ = 0;
  float maxz = -std::numeric_limits<float>::max();
  float minz =  std::numeric_limits<float>::max();

  ForIndex(v, mesh.numVertices()) {
    if (mesh.vertexAt(v)[2] < mesh.vertexAt(min_)[2]) {
      min_ = v;
    }
    maxz = max(maxz, mesh.vertexAt(v)[2]);
    minz = min(minz, mesh.vertexAt(v)[2]);
  }
  float maxh = minz + ((maxz - minz) * layer_thickness / min_thickness);



  // ===================
  // ==== OPTIMIZER ====
  // ===================
  cout << "Creating the model ... ";
  GRBEnv env = GRBEnv();
  AutoPtr<GRBModel> model(new GRBModel(env));
  model->set(GRB_DoubleParam_TimeLimit, compute_time);
  model->set(GRB_IntParam_Threads, nb_threads);
  model->set(GRB_IntParam_Presolve, PRESOLVE);
  model->set(GRB_IntParam_OutputFlag, true);


  // Create variables
  map<uint, GRBVar> h;
  ForIndex(i, mesh.numVertices()) {
    h[i] = model->addVar(0.0, maxh, 0.0, GRB_CONTINUOUS, sprint("z_%03d", i));
  }

/****************************************************************/
/********************     CONSTRAINTS     ***********************/
/****************************************************************/
  if (target_num_layers > -1)
  {
    set<int> surface_vertices;
    ForIndex(t, mesh.numSurfaces()) {
      v3u tri = mesh.surfaceTriangleAt(t);
      ForIndex(i, 3) { surface_vertices.insert(tri[i]); }
    }
    for (auto v : surface_vertices) {
      if (v != i_bottom_shape) {
        model->addConstr(h[v] <= h[i_bottom_shape] + layer_thickness* target_num_layers);
      }
    }
  }
  
  GRBQuadExpr obj = 0.0F;
  ForIndex(t, mesh.numTetrahedrons()) {
    v4u tet = mesh.tetrahedronAt(t);

    if (tetrahedron_volume(mesh,tet) < THRESHOLD_MIN_VOLUME_mm3) {
      continue;
    }

    auto gradient = tetrahedron_gradient(mesh, h, mats, t);
    auto dhdx = gradient[0];
    auto dhdy = gradient[1];
    auto dhdz = gradient[2];


////////////////////////////////////////////////////////////////////
//                                                                //
//                       IGNORE EMPTINESS                         //
//                                                                //
////////////////////////////////////////////////////////////////////
#if IGNORE_EMPTINESS
    if (!mesh.isInside(t)) {
      continue;
    }
#endif

////////////////////////////////////////////////////////////////////
//                                                                //
//                          NO FOLDOVER                           //
//                                                                //
////////////////////////////////////////////////////////////////////
#if CSTRT_FOLDOVER
    if (!mesh.isInside(t) && !fabricable_emptiness)
    {
      model->addConstr(dhdz >= 0.0);
      continue;
    }
#endif

////////////////////////////////////////////////////////////////////
//                                                                //
//                           THICKNESS                            //
//                                                                //
////////////////////////////////////////////////////////////////////
#if CSTRT_THICKNESS
    sl_assert(min_thickness > 0);
    float ratio = max_thickness / min_thickness;
    model->addConstr(dhdz >= 1.0);
    model->addConstr(dhdz <= ratio);
#endif

////////////////////////////////////////////////////////////////////
//                                                                //
//                        COLLISION SLOPE                         //
//                                                                //
////////////////////////////////////////////////////////////////////
#if CSTRT_COLLISION_SLOPE
    double theta = tan(max_theta * M_PI / 180.0);

    model->addConstr(-dhdz <= dhdx / theta);
    model->addConstr( dhdz >= dhdx / theta);
    model->addConstr(-dhdz <= dhdy / theta);
    model->addConstr( dhdz >= dhdy / theta);
#endif


  }

////////////////////////////////////////////////////////////////////
//                                                                //
//                             ANCHOR                             //
//                                                                //
////////////////////////////////////////////////////////////////////
#if CSTRT_ANCHOR_POINT
  model->addConstr(h[min_] == 0);
#elif CSTRT_ANCHOR_LAYER
  // first find out lowest point we would constraint
  ForIndex(i_pts, mesh.numVertices()) {
    if (fabs(mesh.vertexAt(i_pts)[2] - mesh.vertexAt(min_)[2]) < max_thickness / 2.0) {
      model->addConstr(h[i_pts] == (mesh.vertexAt(i_pts)[2] - mesh.vertexAt(min_)[2]));
    }
  }
#endif // CSTRT_ANCHOR_POINT || CSTRT_ANCHOR_LAYER

////////////////////////////////////////////////////////////////////
//                                                                //
//                            FLATTEN                             //
//                                                                //
////////////////////////////////////////////////////////////////////
  set<int> surfaces_to_flatten;
  set<int> relaxed_surfaces;

  ForIndex(s, mesh.numSurfaces()) {
    v3u tri = mesh.surfaceTriangleAt(s);
    if (abs(normals[s][2]) >= normal_threshold) {
      uint a = tri[0], b = tri[1], c = tri[2];
      v3f pa = mesh.vertexAt(a); v3f pb = mesh.vertexAt(b); v3f pc = mesh.vertexAt(c);
      bool bottom = (fabs(pa[2] - bottomz) < 0.05)
                 && (fabs(pb[2] - bottomz) < 0.05)
                 && (fabs(pc[2] - bottomz) < 0.05);
      if (bottom)
      {
        model->addConstr(h[a] - h[b] == 0.0);
        model->addConstr(h[b] - h[c] == 0.0);
      }
      // Ignore flattening for overhanging triangles, unless they are flat already
      else if (normals[s][2] < 0 || normals[s][2] > 0.97)
      { 
        surfaces_to_flatten.insert(s);
      }
    }
  }

/****************************************************************/
/**********************     MAIN LOOP     ***********************/
/****************************************************************/

  AutoPtr<GRBModel> scratch;

  model->update();

  int cnt = 0;
  bool first_pass = true;
  bool last_pass  = false;

  while (!surfaces_to_flatten.empty()) {
    AutoPtr<GRBModel> previous = scratch;

    scratch = AutoPtr<GRBModel>(new GRBModel(*model));

    if (!previous.isNull()) {
      ForIndex(v, mesh.numVertices()) {
        double start = previous->getVarByName(sprint("z_%03d", v)).get(GRB_DoubleAttr_X);
        scratch->getVarByName(sprint("z_%03d", v)).set(GRB_DoubleAttr_Start, start);
      }
    }

    GRBQuadExpr obj = 0.0f;

    {
      // compute total flattened area before filtering
      double tot_flat_area = 0.0;
      for (auto s : surfaces_to_flatten) {
        v3u tri = mesh.surfaceTriangleAt(s);
        tot_flat_area += triangle_area(mesh, tri);
      }

      // filter tiny components
      set<int> filtered_set;
      vector<vector<int>> flattened;
      connected_components(mesh, surfaces_to_flatten, flattened);
      for (const auto& f : flattened) {
        double area = 0.0;
        for (const auto& t : f) {
          v3u tri = mesh.surfaceTriangleAt(t);
          area += triangle_area(mesh, tri);
        }
        if (area > THRESHOLD_AREA_mm2 && f.size() > 1 && (area / tot_flat_area > 0.05)) {
          for (const auto& t : f) {
            filtered_set.insert(t);
          }
        }
      }

      // Remove tiny or isolated areas, and displays the count
      int before = surfaces_to_flatten.size();
      surfaces_to_flatten = filtered_set;
      int after = surfaces_to_flatten.size();
      std::cerr << "removed " << before - after << " tiny or isolated flattened areas." << std::endl;
    }

    double tot_flatten_area = 0.0;
    for (auto s : surfaces_to_flatten) {
      v3u tri = mesh.surfaceTriangleAt(s);
      tot_flatten_area += triangle_area(mesh, tri);
    }

    // ====================
    // ==== flattening ====
    // ====================
    obj_soft_flattening(mesh, tot_flatten_area, h, obj, surfaces_to_flatten);

    // Add alignement objective on each connected component
    double area_check = 0.0;
    bool all_aligned = true;
    if (!previous.isNull()) {
      cerr << Console::white << "alignment objectives" << Console::gray << endl;
      // get components
      vector<vector<int> > flattened;
      connected_components(mesh, surfaces_to_flatten, flattened);
      // sort by area
      vector<pair<double, int> > comp_by_area;
      ForIndex(i,flattened.size()) {
        double comp_area = 0.0;
        for (const auto& t : flattened[i]) {
          v3u tri = mesh.surfaceTriangleAt(t);
          comp_area += triangle_area(mesh, tri);
        }
        comp_by_area.push_back(make_pair(comp_area, i));
      }
      std::sort(comp_by_area.begin(), comp_by_area.end());
      // alignmenent
      double shape_bottom = previous->getVarByName(sprint("z_%03d", i_bottom_shape)).get(GRB_DoubleAttr_X);
      cerr << "shape_bottom = " << shape_bottom << endl;
      // -> foreach component, by decreasing area
      // for (const auto& f : flattened) {
      ForRangeReverse(i,(int)comp_by_area.size()-1,0) {
        const auto& f = flattened[comp_by_area[i].second];
        set<uint> uniquepts;
        double comp_area    = 0.0;
        double comp_area_xy = 0.0;
        bool   comp_flat    = true;
        double comp_non_flatness = 0.0;
        for (const auto& t : f) {
          v3u tri = mesh.surfaceTriangleAt(t);
          comp_area    += triangle_area(mesh, tri);
          comp_area_xy += triangle_area_xy(mesh, tri);
          ForIndex(i, 3) {
            uniquepts.insert(tri[i]);
          }
          // test if flat
          double non_flatness = triangle_non_flatness(tri, previous.raw());
          comp_non_flatness = max(comp_non_flatness, non_flatness);
          if (non_flatness >= thres_flatteness_mm * thres_flatteness_mm) {
            comp_flat = false;
          }
        }

        double nfcomp = component_non_flatness(mesh,f,previous.raw(), max_admissible_non_flatteness_mm);
        comp_flat = (nfcomp < thres_flatteness_mm * thres_flatteness_mm);
        comp_non_flatness = nfcomp;

        area_check += comp_area;

        // compute snapping pos
        double avg = 0.0;

        for (const auto& t : f) {
          v3u tri = mesh.surfaceTriangleAt(t);
          double za = previous->getVarByName(sprint("z_%03d", tri[0])).get(GRB_DoubleAttr_X);
          double zb = previous->getVarByName(sprint("z_%03d", tri[1])).get(GRB_DoubleAttr_X);
          double zc = previous->getVarByName(sprint("z_%03d", tri[2])).get(GRB_DoubleAttr_X);
          avg += triangle_area_xy(mesh, tri) * (za + zb + zc) / 3.0;
        }
        avg /= comp_area_xy;
        avg = avg - shape_bottom;
        double slice = round(avg / layer_thickness) * layer_thickness;
        // align? only if all flat
        bool align   = comp_flat;
        // info
        if (align)     cerr << Console::green; else cerr << Console::white;
        if (comp_flat) cerr << "[flat    ] " << Console::gray;
        else           cerr << Console::yellow << "[non flat] " << Console::gray;
        // cerr << "nonflat: " << comp_non_flatness;
        cerr << " alignement error: " << fabs(avg - slice) << " " << slice <<"-" << avg << " (A = " << comp_area << " w = " << (comp_area / tot_flatten_area) << ")" << endl;
        // do not align => skip
        if (!align) {
          all_aligned = false;
          break; // do not align smaller ones either
        }

        GRBLinExpr lavg = 0.0;
        for (const auto& t : f) {
          v3u tri   = mesh.surfaceTriangleAt(t);
          GRBVar za = scratch->getVarByName(sprint("z_%03d", tri[0]));
          GRBVar zb = scratch->getVarByName(sprint("z_%03d", tri[1]));
          GRBVar zc = scratch->getVarByName(sprint("z_%03d", tri[2]));
          lavg += triangle_area_xy(mesh, tri) * (za + zb + zc) / 3.0;
        }
        lavg /= comp_area_xy;
        GRBLinExpr lexpr = lavg - scratch->getVarByName(sprint("z_%03d", i_bottom_shape)) - slice;

        // weight
        obj += 1.0 * (comp_area / tot_flatten_area) * lexpr*lexpr;


      }
    }
    scratch->update();
    cerr << "area check: expected = " << tot_flatten_area << " found = " << area_check << endl;

    // ====================
    // ==== smoothness ====
    // ====================

#if OBJ_SMOOTHNESS_L2
    obj_smoothness_l2(mesh, h, mats, obj);
#endif // OBJ_SMOOTHNESS_L2

    // ===============
    // ==== slope ====
    // ===============
#if OBJ_SURFACE_SLOPE // normal
    obj_slope(mesh, tot_surface, h, obj, surfaces_to_flatten);
#endif // OBJ_SURFACE_SLOPE


    scratch->setObjective(obj);

    scratch->optimize();


    //std::cout << "Error : " << scratch->getObjectiveError(obj) << std::endl;
    // ==============================================

    {
      TriangleMesh_Ptr tmp(new TriangleMesh_generic<MeshFormat_stl::t_VertexData>(mesh.numVertices(), mesh.numSurfaces()));
      ForIndex(s, mesh.numSurfaces()) {
        tmp->triangleAt(s) = mesh.surfaceTriangleAt(s);
        if (surfaces_to_flatten.find(s) != surfaces_to_flatten.end()) {
          v3u tri = tmp->triangleAt(s);
          std::swap(tri[0], tri[1]);
          tmp->triangleAt(s) = tri;
        }
      }
      ForIndex(v, mesh.numVertices()) {
        tmp->posAt(v) = mesh.vertexAt(v);
        tmp->posAt(v)[2] = scratch->getVarByName(sprint("z_%03d", v)).get(GRB_DoubleAttr_X) + mesh.vertexAt(min_)[2];
      }
      saveTriangleMesh((folder + sprint("/tmp_after_%02d.stl", cnt++)).c_str(), tmp.raw());
    }

    // ==============================================
    if (last_pass) break;
    // ==============================================

    // evaluate
    bool early_stop = true;
    std::vector<std::pair<double, int> > scores;
    std::map<int, double> score_by_tri;
    {
      vector<vector<int> > flattened;
      connected_components(mesh, surfaces_to_flatten, flattened);
      for (const auto& f : flattened) {
#if 1
        {
          double nfcomp = component_non_flatness(mesh, f, scratch.raw(), max_admissible_non_flatteness_mm);
          if (nfcomp < thres_flatteness_mm * thres_flatteness_mm) continue;
        }
#endif
        for (const auto& t : f) {
          v3u   tri = mesh.surfaceTriangleAt(t);
          double non_flatness = triangle_non_flatness(tri, scratch.raw());
          // flat in input?
          if (normals[t][2] > 0.97 || normals[t][2] < -0.97) {
            continue; //  non_flatness = 0.0; // never relaxed
          }
          scores.push_back(make_pair(non_flatness, t));
          score_by_tri[t] = non_flatness;
        }
      }
      std::sort(scores.begin(), scores.end());
    }

    if (!scores.empty()) {
      cerr << Console::cyan;
      cerr << "best  = " << scores.front().first << endl;
      cerr << "worse = " << scores.back().first << endl;
      cerr << "thres = " << thres_flatteness_mm * thres_flatteness_mm << endl;
    }

    cerr << "num flattening 'constraints' " << surfaces_to_flatten.size() << endl;

    // relax some 
    bool none_removed = false;
    bool none_removed_at_all = true;
    /*while (!none_removed)*/ { /* TEST, comment to go back to normal */
      none_removed = true;
      vector<vector<int> > flattened;
      set<int> prev_relaxed = relaxed_surfaces;
      connected_components(mesh, surfaces_to_flatten, flattened);
      for (const auto& f : flattened) { // for each component

#if 1
        {
          double nfcomp = component_non_flatness(mesh, f, scratch.raw(), max_admissible_non_flatteness_mm);
          if (nfcomp < thres_flatteness_mm * thres_flatteness_mm) continue;
        }
#endif

        double worst = 0.0;
        int    worst_t = -1;
        bool   none_removed_comp = true;
        for (const auto& t : f) { // for each triangle
          // flat in input?
          if (normals[t][2] > 0.97 || normals[t][2] < -0.97) {
            continue; // skip
          }
          // track worst
          if (score_by_tri[t] > worst) {
            worst = score_by_tri[t];
            worst_t = t;
          }
          // on border?
          bool onborder = false;
          vector<uint> neighs = mesh.tri_neighbours(t);
          ForIndex(n, neighs.size()) {
            // check if neighbor has been relaxed
            if (prev_relaxed.find(neighs[n]) != prev_relaxed.end()) {
              // there is at least one relaxed neighbor nearby
              onborder = true;
              break;
            }
          }
          if (onborder && score_by_tri[t] >= thres_flatteness_mm * thres_flatteness_mm) {
            // kill 
            none_removed = false;
            none_removed_comp = false;
            surfaces_to_flatten.erase(t);
            relaxed_surfaces.insert(t);
          }
        }
        if (none_removed_comp && worst_t > -1) {
          if (score_by_tri[worst_t] >= thres_flatteness_mm * thres_flatteness_mm) {
            // kill
            none_removed = false;
            surfaces_to_flatten.erase(worst_t);
            relaxed_surfaces.insert(worst_t);
          }
        }
      }
      if (!none_removed) {
        none_removed_at_all = false;
      }
    }

    if (!scores.empty()) {

      // safety: kill at least the worst (unless it is below theshold)
      if (none_removed_at_all && scores.back().first > thres_flatteness_mm*thres_flatteness_mm) {
        cerr << Console::yellow << "WARNING: relaxation safety triggered" << Console::gray << endl;
        //sl_assert(false); // SL: no longer supposed to happen
        none_removed_at_all = false;
        surfaces_to_flatten.erase(scores.back().second);
        relaxed_surfaces.insert(scores.back().second);
      }

      if (!first_pass) {

        if (scores.back().first < thres_flatteness_mm*thres_flatteness_mm && none_removed_at_all) {
          // termination on success!
          cerr << Console::gray;
          last_pass = true && all_aligned;
        }

      }
    } else {
      last_pass = true;
    }
    first_pass = false;

    cerr << "num flattening 'constraints' after relaxation " << surfaces_to_flatten.size() << endl;
    cerr << Console::gray;

  }

  model = scratch;

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  ////////////   ALIGNMENT CHECK  ///////////////
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////

  vector<vector<int> > flattened;
  connected_components(mesh, surfaces_to_flatten, flattened);
  double shape_bottom = model->getVarByName(sprint("z_%03d", i_bottom_shape)).get(GRB_DoubleAttr_X);
  for (const auto& f : flattened) {
    set<uint> uniquepts;
    double area = 0.0;
    for (const auto& t : f) {
      v3u tri = mesh.surfaceTriangleAt(t);
      ForIndex(i, 3) {
        uniquepts.insert(tri[i]);
      }
      area += triangle_area(mesh, tri);
    }
    double avg = 0.0;
    for (auto v : uniquepts) {
      double vh = model->getVarByName(sprint("z_%03d", v)).get(GRB_DoubleAttr_X);
      avg += vh;
    }
    avg /= (double)uniquepts.size();
    double slice = shape_bottom + round((avg - shape_bottom) / layer_thickness)*layer_thickness;
    cerr << "alignement error: " << fabs(avg - slice) << "(area: " << area << ")" << endl;
  }

  //cerr << "=============== POSTPROCESS =============== " << endl;
  //
  //AutoPtr<SLRModel<double>> post = postprocess(mesh,env,model,maxh,i_bottom_shape,surfaces_to_flatten);
  AutoPtr<GRBModel> post = model;

  // ==============================================

  // get the result
  ForArray(h, i) {
    GRBVar vh = post->getVarByName(sprint("z_%03d", i));
    n_pts[i] = float(vh.get(GRB_DoubleAttr_X) + mesh.vertexAt(min_)[2]);
  }

  // ==============================================

  map<uint, float> z;
  for (uint i = 0; i < n_pts.size(); i++) {
    z[i] = n_pts[i];
  }

  double lmin = std::numeric_limits<double>::max();
  double lmax = -std::numeric_limits<double>::max();

  ForIndex(i_srf, mesh.numSurfaces()) {
    ForIndex(i_ver, 3) {
      double z = n_pts[mesh.surfaceTriangleAt(i_srf)[i_ver]];
      lmin = std::min(lmin, z);
      lmax = std::max(lmax, z);
    }
  }
  cerr << "Vertical extent: " << lmax - lmin << endl;

  // ==============================================

  int plot_nb_cstr = surfaces_to_flatten.size();
  int nb_slices = floor((lmax - lmin) / layer_thickness);
  double plot_slope_obj = 0.0;
  obj_slope(mesh, tot_surface, z, plot_slope_obj, surfaces_to_flatten);
  double plot_slope_obj_incl_flat = 0.0;
  auto tmpVar = set<int>();
  obj_slope(mesh, tot_surface, z, plot_slope_obj_incl_flat, tmpVar);

  std::ofstream out;
  out.open(folder + "/objective_" + filename + ".csv", std::ios::app);
  out << normal_threshold << "," << obj_angle << "," << nb_slices << "," << time.elapsed() << "," << plot_slope_obj << "," << plot_slope_obj_incl_flat << endl;

  // ==============================================

  return model->get(GRB_IntAttr_SolCount);
}

// --------------------------------------------------------------

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
  string path_models = SRC_PATH;
  max_theta = 30.0;
  obj_angle = 0.0;
  compute_time = 6*600; // Default is 60 min

  // command line
  TCLAP::CmdLine cmd("", ' ', "1.0");
  TCLAP::UnlabeledValueArg<std::string>  dirArg("d", "dir" , true, "models\\", "directory name");
  TCLAP::UnlabeledValueArg<std::string> fileArg("f", "file", true, "filename", "file name"     );

  TCLAP::ValueArg<float>       tauArg("t", "tau"             , "layer thickness (mm)"            , false,                    0.0f, "float");
  TCLAP::ValueArg<float>     thetaArg("" , "theta"           , "Maximum printing angle (degrees)", false,                   30.0f, "float");
  TCLAP::ValueArg<float> thresholdArg("n", "normal-threshold", "Normal threshold"                , false,                    0.0f, "float");
  TCLAP::ValueArg<float>   layerThArg("l", "layer-thickness" , "Layer thickness (mm)"            , false, default_layer_thickness, "float");
  TCLAP::ValueArg<float>       phiArg("" , "phi"             , "deprecated"                      , false,                     0.0, "float");
  TCLAP::ValueArg<float>     ratioArg("r", "ratio"           , "ratio"                           , false,                     6.0, "float");

  TCLAP::ValueArg<int>    numLayerArg("" , "numlayers"   , "Target number of layers"                               , false,  -1, "int");
  TCLAP::ValueArg<int>      threadArg("T", "threads"     , "number of threads used for optimisation, default: 4"   , false,   4, "int");
  TCLAP::ValueArg<int>        timeArg("c", "compute-time", "max compute time (s) for each iteration, default: 600s", false, 600, "int");
  
  TCLAP::SwitchArg        fabemptyArg("" , "fabempty"    , "fabricable emptyness"                                  , true);
  

  // file and folder
  cmd.add(dirArg);
  cmd.add(fileArg);

  // angles and thresholds
  cmd.add(tauArg);
  cmd.add(thetaArg);
  cmd.add(phiArg);

  // thicknesses and ratio
  cmd.add(layerThArg);
  cmd.add(numLayerArg);
  cmd.add(ratioArg);

  cmd.add(thresholdArg);
  cmd.add(fabemptyArg);
  
  // computation
  cmd.add(timeArg);
  cmd.add(threadArg);

  // Parsing
  cmd.parse(argc, argv);

  path_models          = dirArg.getValue();
  filename             = fileArg.getValue();
  
  if ( thetaArg.isSet()) max_theta = thetaArg.getValue();

  layer_thickness   = layerThArg.getValue();
  target_num_layers = numLayerArg.getValue();

  max_thickness = layer_thickness;
  min_thickness = max_thickness / ratioArg.getValue();

  normal_threshold = thresholdArg.isSet()
                     ? thresholdArg.getValue()
                     : cos(max_theta * M_PI / 180);
  fabricable_emptiness = fabemptyArg.getValue();

  if (  timeArg.isSet()) compute_time = (float)timeArg.getValue();
  if (threadArg.isSet())   nb_threads = threadArg.getValue();

  /////////////////////////////

  folder = path_models + removeExtensionFromFileName(filename);

  try {

    TetMesh* mesh(TetMesh::load((path_models + filename).c_str()));
    cerr << Console::white << "Mesh has " << mesh->numTetrahedrons() << " tets." << Console::gray << std::endl;

    if (numLayerArg.isSet()) {
      float bottomz = std::numeric_limits<float>::max();
      int   i_bottom_shape = -1;
      ForIndex(t, mesh->numSurfaces()) {
        v3u tri = mesh->surfaceTriangleAt(t);
        v3f p[3] = { mesh->vertexAt(tri[0]), mesh->vertexAt(tri[1]), mesh->vertexAt(tri[2]) };
        ForIndex(i, 3) {
          if (p[i][2] < bottomz) {
            bottomz = p[i][2];
          }
        }
      }

      bool toMuch = bottomz > (layer_thickness * max_thickness / min_thickness) * target_num_layers;
      bool toFew = layer_thickness * target_num_layers < bottomz;
      if (toMuch || toFew) {
        cerr << Console::red << "Targeted number of layers unreachable, ";
        if (toMuch)
          cerr << "to much, " << mesh->getBBox().extent()[2] << " > " << (layer_thickness * max_thickness / min_thickness) * target_num_layers;
        if (toFew)
          cerr << "to few, " << mesh->getBBox().extent()[2] << " < " << layer_thickness * target_num_layers;

        cerr << Console::gray << std::endl;
        return 0;
      }
    }

    {
      // load triangle mesh
      createDirectory(folder.c_str());
      mesh->save((folder + "/before.stl").c_str(), mesh);

      Array<Tuple<double, 9>> mats;
      {
        // prepare tetrahedrons
        cout << "Compute all matrices... ";
        cerr << mesh->numTetrahedrons() << endl;
        mats.allocate(mesh->numTetrahedrons());
        ForIndex(t, mesh->numTetrahedrons()) {
          mats[t] = mesh->getGradientMatrix(t);
        }
        saveArray(mats, (folder + "/tetmats").c_str());
        cout << "done !" << endl;
      }

      // optimize
      Array<float> displ(mesh->numVertices());



      try {
        displ.fill(0.0f);
        gurobi_opt(*mesh, mats, displ);
      } catch (const GRBException& e) {
        cerr << Console::red << "Optimizer error " << e.getErrorCode() << ": " << e.getMessage() << Console::gray << endl;
        return -1;
      }

      cout << Console::green << "Solution found !" << endl;

      TetMesh tmp(mesh);
      ForArray(displ, i) {
        tmp.vertexAt(i)[2] = displ[i];
      }

      tmp.save((folder + "/after.stl").c_str(), &tmp);
      saveArray(displ, (folder + "/displacements").c_str());

      ///////////////////////////////////////

      if (thresholdArg.isSet() || layerThArg.isSet()) {

        std::ostringstream angle;
        angle.precision(0);
        angle << std::fixed << obj_angle;

        std::ostringstream cstr_angle;
        cstr_angle.precision(0);
        cstr_angle << std::fixed << max_theta;

        std::ostringstream normal;
        normal.precision(3);
        normal << std::fixed << normal_threshold;

        std::ostringstream lthick;
        lthick.precision(2);
        lthick << std::fixed << layer_thickness;

        std::string postfix = angle.str() + "_" + cstr_angle.str() + "_" + normal.str() + "_" + lthick.str();

        saveArray(displ, (folder + "/displacements_" + postfix).c_str());
        tmp.save((folder + "/after_" + postfix + ".stl").c_str(), &tmp);
      }

    }
  } catch (const GRBException& e) {
    cerr << Console::red << "Optimizer error " << e.getErrorCode() << ": " << e.getMessage() << Console::gray << endl;
  } catch (Fatal& e) {
    cerr << "[ERROR] " << e.message() << endl;
  }
  cerr << Console::gray;
}

// --------------------------------------------------------------
