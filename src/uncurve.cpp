/*
This work and all associated files are under the

     GNU AFFERO GENERAL PUBLIC LICENSE
        Version 3, 19 November 2007

A copy of the license full text is included in
the distribution, please refer to it for details.

(c) Jimmy Etienne and Sylvain Lefebvre
*/

#include <algorithm>
#include <LibSL/LibSL.h>
#include <LibSL/LibSL.h>
#include <LibSL/Memory/ArrayTools.h>
#include <tclap/CmdLine.h>

#include <iostream>
#include <queue>
#include <iomanip>
#include <set>
#include <iterator>
#include <unordered_map>

#include "Model.h"

#include "gcode.h"
#include "TetMesh.h"
#include "MeshFormat_msh.h"

#include "thicknesses.h"

// --------------------------------------------------------------
using namespace std;

const float retract_e_length_mm = 6.0f;
const float min_speed_mm_sec    = 10.0f;
const float max_speed_mm_sec    = 30.0f;
const double e_dampening        = 0.025;
const float c_delta_centering = 75.0f;

bool um2       = false;
bool a8        = false;
bool deltaP    = false;
bool do_reflow = true;
bool do_sim    = false;

// --------------------------------------------------------------

// Commandline parameters

string filename;

float layer_thickness = default_layer_thickness;

// --------------------------------------------------------------

#define ScTP(a, b, c) dot(a, cross(b,c))

// --------------------------------------------------------------

float tetrahedron_volume(const Tuple<v3f, 4>& tet)
{

  v3f vab = tet[1] - tet[0];
  v3f vac = tet[2] - tet[0];
  v3f vad = tet[3] - tet[0];

  // ScTP computes the scalar triple product
  float v = 1.f / 6.f * ScTP(vab, vac, vad);

  return abs(v);
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

  // ScTP computes the scalar triple product
  float va = 1.f / 6.f * ScTP(vbp, vbd, vbc);
  float vb = 1.f / 6.f * ScTP(vap, vac, vad);
  float vc = 1.f / 6.f * ScTP(vap, vad, vab);
  float vd = 1.f / 6.f * ScTP(vap, vab, vac);
  float v  = 1.f / 6.f * ScTP(vab, vac, vad);

  float v_ = abs(v);
  v4f coeffs = v4f(va / v_, vb / v_, vc / v_, vd / v_);
  return coeffs;
}

// --------------------------------------------------------------

float triangle_area(const Tuple<v3f, 3>& tri)
{

  v3f vab = tri[1] - tri[0];
  v3f vac = tri[2] - tri[0];

  float v = 1.f / 2.f * length(cross(vab, vac));

  return abs(v);
}

// --------------------------------------------------------------

float triangle_area(TetMesh& mesh, const v3u& tri)
{
  return triangle_area(Tuple<v3f, 3>(mesh.vertexAt(tri[0]), mesh.vertexAt(tri[1]), mesh.vertexAt(tri[2])));
}

// --------------------------------------------------------------

bool sameSide(const v3f& v1, const v3f& v2, const v3f& v3, const v3f& v4, const v3f& pt)
{
  v3f  normal = cross(v2 - v1, v3 - v1);
  float dotV4 = dot(normal, v4 - v1);
  float dotPt = dot(normal, pt - v1);
  return (sign(dotV4) * sign(dotPt)) >= -1e-6f;
}

// --------------------------------------------------------------

inline bool isInTetrahedron(const Tuple<v3f, 4>& tet, const v3f& pt)
{
  float vol = tetrahedron_volume(tet);
  if (vol < 0.01f) { // see curvislice thres_min_tet_volume_mm3
    return false;
  }
  float eps = 1e-9f;
  v4f coefs = tetrahedron_barycenter_coefs(tet, pt);
  for (int i = 0; i < 4; ++i) {
    if (coefs[i] < -eps || coefs[i] > 1.f + eps) {
      return false;
    }
  }
  return true;

  /*
  return sameSide(tet[0], tet[1], tet[2], tet[3], pt) &&
         sameSide(tet[1], tet[2], tet[3], tet[0], pt) &&
         sameSide(tet[2], tet[3], tet[0], tet[1], pt) &&
         sameSide(tet[3], tet[0], tet[1], tet[2], pt);
    */
}

// --------------------------------------------------------------

template<typename T_Type>
T_Type tetrahedral_interpolation(const Tuple<T_Type, 4>& tet, v4f& coefs)
{

  T_Type vab = tet[1] - tet[0];
  T_Type vac = tet[2] - tet[0];
  T_Type vad = tet[3] - tet[0];

  T_Type interp = coefs[0] * tet[0]
    + coefs[1] * tet[1]
    + coefs[2] * tet[2]
    + coefs[3] * tet[3];

  return interp;
}

// --------------------------------------------------------------

bool                 init = false;
Array3D<vector<int>> tet_location;
v3f                  split; // bad name, but just for here
AABox                meshbox;

AABox getBBox(Tuple<v3f, 4> tet)
{
  AABox bbox;
  ForIndex(i, 4) {
    bbox.addPoint(tet[i]);
  }
  return bbox;
}

int find_containing_tet(TetMesh* mesh, v3f& pt, bool use_in_and_out)
{
  int res = 100;
  if (!init) {
    init = true;

    meshbox = mesh->getBBox();

    tet_location.allocate(res, res, res);
    tet_location.fill(vector<int>());

    split = meshbox.maxCorner() - meshbox.minCorner();
    split[0] /= res; split[1] /= res; split[2] /= res;

    ForIndex(itet, mesh->numTetrahedrons()) {

      Tuple<v3f, 4> tet;
      ForIndex(i, 4) {
        tet[i] = mesh->vertexAt(mesh->tetrahedronAt(itet)[i]);
      }

      AABox tetbox = getBBox(tet);

      v3i start = v3i(floor((tetbox.minCorner() - meshbox.minCorner()) / split));
      v3i end = v3i(ceil((tetbox.maxCorner() - meshbox.minCorner()) / split));

      for (int x = start[0]; x <= end[0]; ++x) {
        for (int y = start[1]; y <= end[1]; ++y) {
          for (int z = start[2]; z <= end[2]; ++z) {
            tet_location.at<Clamp>(x, y, z).emplace_back(itet);
          }
        }
      }
    }
  }

  v3i pos = v3i(round((pt - meshbox.minCorner()) / split));
  ForIndex(i, 3) pos[i] = min(max(pos[i], 0), res - 1);

  const vector<int>& tets = tet_location.at<Clamp>(pos[0], pos[1], pos[2]);

  // traverses the entire list, return the first solid, or the first encountered empty

  int itet_empty_containing = -1;

  for (int itet : tets) {

    Tuple<v3f, 4> tet;

    ForIndex(i, 4) {
      tet[i] = mesh->vertexAt(mesh->tetrahedronAt(itet)[i]);
    }

    if (isInTetrahedron(tet, pt)) {
      if (mesh->isInside(itet)) {
        return itet;
      } else {
        itet_empty_containing = itet;
      }
    }
  }

  return itet_empty_containing;
}

// --------------------------------------------------------------

enum e_step_type { e_Prime, e_Retract, e_Print, e_Travel, e_Ironing };

typedef struct {
  e_step_type type;
  v3f         xyz;
  float       e;
  float       f;
} t_step;

#define SAMPLING 0.8f

// --------------------------------------------------------------


void reflow(std::vector<t_step>& _steps)
{
  const double max_e_delta = 1.0;

  // optimize flow
  std::vector<SLRVar<double>> vars;
  vars.resize(2*_steps.size());

  // configure gurobi env
  AutoPtr<SLRModel<double>> model(new SLRModel<double>());

  model->setTimeLimit(600);
  //model->setNbThread(8);
  model->setPresolve(-1);
  model->printDebug(true);

  double e_start = _steps.front().e;

  // vars
  ForIndex(i,_steps.size()) {
    if (_steps[i].type == e_Print) {
      vars[i]               = model->addVar(0.0, (_steps.back().e - e_start) * 2.0f, 0.0, sprint("e_%03d", i));
      vars[i+_steps.size()] = model->addVar(0.0, (_steps.back().e - e_start) * 2.0f, 0.0, sprint("a_%03d", i));
    }
  }

  // constraints
#if 1
  SLRVar<double> e_i_m_1;
  SLRVar<double> a_i_m_1;
  bool first = true;
  v3f prev_pos;
  ForIndex(i, _steps.size()) {
    if (_steps[i].type == e_Print) {
      auto e_i = vars[i];
      auto a_i = vars[i + _steps.size()];
      if (!first) {
        // limit gradient on E
        model->addConstr(e_i <= e_i_m_1 + 1.0f * SAMPLING       );
        model->addConstr(       e_i_m_1 - 1.0f * SAMPLING <= e_i);
        // compute a_i
        float w = min(1.0f, length(_steps[i].xyz - prev_pos) / SAMPLING);
        model->addConstr(a_i == a_i_m_1 + w * e_dampening * (e_i - a_i_m_1));
      } else {
        model->addConstr(vars[i] == 0.0); // e_0
        model->addConstr(vars[i+_steps.size()] == 0.0); // a_0
        first = false;
      }
      prev_pos = _steps[i].xyz;
      e_i_m_1 = e_i;
      a_i_m_1 = a_i;
    }
  }
#endif

  // objective for reflow
  double      w_min = 0.0001;
  SLRExpr<double> obj   = 0.0f;
  unordered_map<int, double> weights;
  Console::progressTextInit((int)_steps.size());
  ForIndex(i, _steps.size()) {
    Console::progressTextUpdate(i);
    if (_steps[i].type == e_Print) {
      auto   a_i = vars[i + _steps.size()];
      double E_i = (_steps[i].e - e_start);
      obj += (a_i - E_i) * (a_i - E_i);
    }
  }
  Console::progressTextEnd();

  LIBSL_TRACE;

  model->setObjective(obj);

  LIBSL_TRACE;

  model->optimize();

  // set result
  float e = 0.0f;
  ForIndex(i, _steps.size()) {
    if (_steps[i].type == e_Print) {
      float v = (float)vars[i].get();
      e       = v; // _steps[i].e;
    }
    _steps[i].e = e;
  }

}


// --------------------------------------------------------------

void inv_curv(TetMesh* mesh, string filepath)
{

  // inverse displacement
  Array<float> h;
  loadArray(h, (filepath + "/displacements").c_str());

  TetMesh* mesh_crv = new TetMesh(mesh);
  ForArray(h, i) {
    mesh_crv->vertexAt(i)[2] = h[i];
  }

  Array<Tuple<double, 9>> mats;
  loadArray(mats, (filepath + "/tetmats").c_str());

  std::string gcode = loadFileIntoString((filepath + ".gcode").c_str());

  gcode_start(gcode.c_str());

  if (fabs(gcode_layer_thickness() - layer_thickness) > 1e-3f) {
    cerr << Console::red;
    cerr << "GCode thickness : " << gcode_layer_thickness() << endl;
    cerr << "Layer thickness : " << layer_thickness << endl;
    cerr << "=====================================================\n";
    cerr << "=====================================================\n";
    cerr << " GCode layer thickness does not match compiled parameter\n";
    cerr << "=====================================================\n";
    cerr << "=====================================================\n";
    cerr << Console::gray;

    layer_thickness = gcode_layer_thickness();
    exit(-1);
  }

  v4f offset(gcode_offset(), 0.0f); // this brings the gcode 'z=0' to object 'h=0'

  float zmin_object = std::numeric_limits<float>::max();
  float hmin_object = std::numeric_limits<float>::max();
  ForIndex(t, mesh->numTetrahedrons()) {
    if (mesh->isInside(t)) {
      ForIndex(i, 4) {
        zmin_object = min(zmin_object, mesh->vertexAt(mesh->tetrahedronAt(t)[i])[2]);
        hmin_object = min(hmin_object, h[mesh->tetrahedronAt(t)[i]]);
      }
    }
  }

  v2f centering = v2f(0.0f);
  if (deltaP) {
    centering = - 75.0f;
  }

  ////////////////////////////

  v4f pos, prev_pos;

  gcode_advance();
  prev_pos = gcode_next_pos() + offset;

  std::vector< t_step > steps;

  float min_thick = std::numeric_limits<float>::max();
  float max_thick = -std::numeric_limits<float>::max();

  float z_clamp = -std::numeric_limits<float>::max();
  float z_last = -std::numeric_limits<float>::max();

  bool  traveling = false;

  float total_e = 0.0f;

  // gather steps
  while (gcode_advance()) {

    if (gcode_do_prime()) {

      steps.push_back(t_step());
      steps.back().e = total_e;
      steps.back().f = 0.0f;
      steps.back().xyz = v3f(0.0f,0.0f,z_last);
      steps.back().type = e_Prime;

      traveling = false;
      continue;

    } else if (gcode_do_retract()) {

      steps.push_back(t_step());
      steps.back().e = total_e;
      steps.back().f = 0.0f;
      steps.back().xyz = v3f(0.0f, 0.0f, z_last);
      steps.back().type = e_Retract;

      traveling = true;
      continue;
    }

    // read gcode
    pos = gcode_next_pos() + offset;
    float f = gcode_speed();

    // get distance
    v3f  delta   = v3f(pos - prev_pos);
    float dist   = length(delta);
    int nb_steps = max(2, int(dist / SAMPLING));

    ForIndex(step, nb_steps) {

      float a      = step / float(nb_steps);
      float next_a = (step + 1) / float(nb_steps);

      // field sampling point ; note the half layer thickness shift to sample on the slicing plane
      v4f A = (a * pos + (1.f - a) * prev_pos) - v4f(0, 0, layer_thickness / 2.0f, 0.0f);
      v4f B = (next_a * pos + (1.f - next_a) * prev_pos) - v4f(0, 0, layer_thickness / 2.0f, 0.0f);

      v4u itet;
      Tuple<v3f, 4> tet, tet_crv;
      v4f coefs;

      bool A_inside = false;

      int i;
      v3f vectorA(A);
      if ((i = find_containing_tet(mesh_crv, vectorA, true)) != -1) {
        itet = mesh->tetrahedronAt(i);
        for (int i = 0; i < 4; i++) { tet[i] = mesh->vertexAt(itet[i]); tet_crv[i] = mesh_crv->vertexAt(itet[i]); }
        coefs = tetrahedron_barycenter_coefs(tet_crv, v3f(A));
        A[2] = tetrahedral_interpolation(tet, coefs)[2];

        A_inside = (mesh->isInside(i));

      } else {
        // if no tet found, do nothing
        continue;
      }

      float ha = h[itet[0]];
      float hb = h[itet[1]];
      float hc = h[itet[2]];
      float hd = h[itet[3]];

      float grad_z_a = (ha - hd) * (float)mats[i][6]
        + (hb - hd) * (float)mats[i][7]
        + (hc - hd) * (float)mats[i][8];

      if (A_inside) {
        min_thick = min(min_thick, layer_thickness / grad_z_a);
        max_thick = max(max_thick, layer_thickness / grad_z_a);
      }

      A[2] += layer_thickness * 0.5f / grad_z_a;

      v3f vectorB(B);
      if ((i = find_containing_tet(mesh_crv, vectorB, true)) != -1) {
        itet = mesh->tetrahedronAt(i);
        for (int i = 0; i < 4; i++) { tet[i] = mesh->vertexAt(itet[i]); tet_crv[i] = mesh_crv->vertexAt(itet[i]); }
        coefs = tetrahedron_barycenter_coefs(tet_crv, v3f(B));
        B[2] = tetrahedral_interpolation(tet, coefs)[2];
      } else {
        // if no tet found, do nothing
       // sl_assert(false);
        continue;
      }

      ha = h[itet[0]];
      hb = h[itet[1]];
      hc = h[itet[2]];
      hd = h[itet[3]];

      float grad_z_b = (ha - hd) * (float)mats[i][6]
        + (hb - hd) * (float)mats[i][7]
        + (hc - hd) * (float)mats[i][8];

      B[2] += layer_thickness * 0.5f / grad_z_b;

      float delta_E   = B[3] - A[3];
      float delta_pos = length(v2f(B) - v2f(A));

      float grad = 0.5f * grad_z_a + 0.5f * grad_z_b;
      // float grad = max( grad_z_a , grad_z_b );

      float new_e = delta_E / grad;
      total_e   += new_e;

      float ratio = delta_E / new_e;

      if (!um2) {
        f = gcode_speed() * grad;
        f = min(max_speed_mm_sec, max(min_speed_mm_sec, f)) * 60.f;
      } else {
        f = min_speed_mm_sec * 60.0f;
      }

      B[0] -= offset[0];
      B[1] -= offset[1];
      // adjust based on offset input/deformed object
      B[2] += -zmin_object;

      float z = B[2];
      z_last = z;
      if (traveling) {
        sl_assert(!gcode_is_ironing());
        z = max(B[2], z_clamp);
      }

      steps.push_back(t_step());
      steps.back().e = total_e;
      steps.back().f = f;
      steps.back().xyz = v3f(B[0],B[1], max(0.0f, z));
      steps.back().type = traveling ? e_Travel : (gcode_is_ironing() ? e_Ironing : e_Print);
    }

    // next
    prev_pos = pos;

  }
  cerr << "num steps: " << steps.size() << endl;

  ////////////////////////////////////////////////////////

  if (do_reflow) {
#if 1

    reflow(steps);

#else
    std::vector< t_step > new_steps;

    // split into paths and solve each
    double e_prev = 0.0;
    int i = 0;
    while (i < steps.size()) {

      std::vector< t_step > sub_steps;
      while (steps[i].type != e_Retract && i < steps.size()) {
        sub_steps.push_back(steps[i]);
        i++;
      }
      sub_steps.push_back(steps[i]);
      i++;

      reflow(sub_steps);

      for (auto& step : sub_steps) {
        step.e += e_prev;
      }
      e_prev = sub_steps.back().e;

      new_steps.insert(new_steps.end(), sub_steps.begin(), sub_steps.end());
    }

    steps = new_steps;
#endif
  }

  ////////////////////////////////////////////////////////

  std::ofstream out_gcode((filepath + ".gcode").c_str(), ios::out | ios::trunc);
  if (!out_gcode) {
    std::cerr << "Erreur � l'ouverture !" << endl;
    return;
  }

  float actual_e = 0.0f;
  prev_pos = v4f(steps.front().xyz);
  ForIndex(s,steps.size()) {

    const auto& step = steps[s];

    if (step.type == e_Prime) {
      if (um2) {
        out_gcode << "G1 Z" << std::setprecision(5) << step.xyz[2] << " ; prime (z correction)" << std::endl;
        out_gcode << "G11" << std::endl;
      } else {
        out_gcode << "G1 E" << std::setprecision(7) << step.e << " ; prime" << std::endl;
      }
      continue;
    } else if (step.type == e_Retract) {
      if (um2) {
        out_gcode << "G10" << std::endl;
      } else {
        out_gcode << "G1 E" << std::setprecision(7) << step.e - retract_e_length_mm << " ; retract" << std::endl;
      }
      continue;
    } else {

      if (step.type == e_Print) {

        float w = min(1.0f, length(step.xyz - v3f(prev_pos)) / SAMPLING);
        prev_pos = v4f(step.xyz);
        actual_e = actual_e + (step.e - actual_e) * (float)e_dampening * w;

        out_gcode << "G1 X" << step.xyz[0] + centering[0] << " Y" << step.xyz[1] + centering[1] << " Z" << max(0.0f, step.xyz[2]);
        if (do_sim) {
          out_gcode << " E" << std::setprecision(7) << actual_e;
        } else {
          out_gcode << " E" << std::setprecision(7) << step.e;
        }
        out_gcode << " F" << std::setprecision(5) << step.f << endl;

      } else {

        // lift up if ironing
        float z_lift = 0.0;
        if (step.type == e_Ironing) {
          const float nozzle_padding = 1.0f;
          if (s > 0 && s + 1 < steps.size()) { // if enough points around
            v3f m1 = steps[s - 1].xyz;
            v3f p1 = steps[s + 1].xyz;
            v3f delta = p1 - m1;
            if (sqLength(delta) > 1e-6f) {
              float z_delta = delta[2];
              float agl = asin(z_delta / length(delta));
              z_lift = abs(tan(agl) * nozzle_padding);
            }
          }
        }

        out_gcode << "G0 X" << step.xyz[0] + centering[0] << " Y" << step.xyz[1] + centering[1] << " Z" << max(0.0f, step.xyz[2] + z_lift)
            << " F" << std::setprecision(5) << step.f << endl;
      }
    }
  }

  out_gcode.close();

  cerr << "min/max thickness encountered: " << min_thick << "/" << max_thick << endl;
}

// --------------------------------------------------------------

void inv_deform(TetMesh* mesh, string path)
{
  // inverse displacement
  Array<float> h;
  loadArray(h, (path + "/displacements").c_str());

  TetMesh* mesh_crv = new TetMesh(mesh);
  ForArray(h, i) {
    mesh_crv->vertexAt(i)[2] = h[i];
  }

  AABox bbox = mesh_crv->getBBox();

  TriangleMesh *stl = loadTriangleMesh((path + ".stl").c_str());

  ForIndex(v, stl->numVertices()) {

    MeshFormat_stl::t_VertexData *d = (MeshFormat_stl::t_VertexData*)stl->vertexDataAt(v);

    v4u itet;
    Tuple<v3f, 4> tet, tet_crv;
    v4f coefs;

    int i;
    if ((i = find_containing_tet(mesh_crv, d->pos, true)) != -1) {
      itet = mesh->tetrahedronAt(i);
      for (int i = 0; i < 4; i++) { tet[i] = mesh->vertexAt(itet[i]); tet_crv[i] = mesh_crv->vertexAt(itet[i]); }
      coefs = tetrahedron_barycenter_coefs(tet_crv, d->pos);
      d->pos[2] = tetrahedral_interpolation(tet, coefs)[2];;
    } else {
      if (!bbox.contain(d->pos)) {
        cerr << 'o';
      }
    }
  }

  MeshFormat_stl meshformat;
  meshformat.save((path + "/uncurved.stl").c_str(), stl);
}

// --------------------------------------------------------------

int main(int argc, char **argv)
{
  string path = "";

  // command line
  TCLAP::CmdLine cmd("", ' ', "1.0");
  TCLAP::UnlabeledValueArg<std::string>  dirArg("d", "path", true, ".", "path");

  TCLAP::ValueArg<float>   layerThArg("l", "layer-thickness", "Layer thickness (mm)", false, default_layer_thickness, "float");

  TCLAP::SwitchArg invArg("g", "gcode", "generate curved gcode", true);
  TCLAP::SwitchArg stlArg("s", "stl", "generate stl from msh with deformations", true);
  TCLAP::SwitchArg um2Arg("", "um2", "generate GCode for UM2", false);
  TCLAP::SwitchArg a8Arg("", "a8", "generate GCode for A8", false);
  TCLAP::SwitchArg deltaArg("", "delta", "generate GCode for the delta printer", false);
  TCLAP::SwitchArg rflwArg("", "reflow", "apply reflow", false);
  TCLAP::SwitchArg simArg("", "sim", "'simulate' nozzle pressure", false);

  cmd.add(invArg);
  cmd.add(stlArg);
  cmd.add(um2Arg);
  cmd.add(a8Arg);
  cmd.add(deltaArg);
  cmd.add(rflwArg);
  cmd.add(simArg);

  cmd.add(layerThArg);
  cmd.add(dirArg);

  cmd.parse(argc, argv);

  path          = dirArg  .getValue();
  max_thickness = layerThArg.getValue();
  um2           = um2Arg  .getValue();
  a8            = a8Arg   .getValue();
  deltaP        = deltaArg.getValue();
  do_reflow     = rflwArg .getValue();
  do_sim        = simArg  .getValue();

  /////////////////////////////

  layer_thickness = max_thickness;
  min_thickness = max_thickness / 6.0f;

  path = removeExtensionFromFileName(path);

  try {

    TetMesh* mesh(TetMesh::load((path + ".msh").c_str()));

    if (invArg.getValue() == false) {
      inv_curv(mesh, path);
    } else if (stlArg.getValue() == false) {
      inv_deform(mesh, path);
    } else {
      std::cerr << "nothing to do ..." << std::endl;
    }
  } catch (SLRException& e) {
    cerr << "[ERROR] " << e.getMessage() << endl;
  } catch (Fatal& e) {
    cerr << "[ERROR] " << e.message() << endl;
  }
}

// --------------------------------------------------------------
