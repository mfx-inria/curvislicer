#include "TetViz.h"
#include "shade.h"
#include "slicerror.h"

#include <iostream>
#include <ctime>
#include <cmath>

#include <LibSL/Mesh/Mesh.h>
#include <LibSL/Mesh/MeshFormat_msh.h>

#include <LibSL/LibSL.h>
#include <LibSL/LibSL_gl4.h>
#include <LibSL/GPUHelpers/GPUHelpers.h>

using namespace LibSL::Mesh;
using namespace LibSL;
// --------------------------------------------------------------

using namespace std;

// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------

double thickness = 0.6; ///////// SL:WARNING!!!! hard coded thickness (also elsewehre)

// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------

LIBSL_WIN32_FIX

GLUX_LOAD(GL_VERSION_2_0)

// --------------------------------------------------------------

#include <tclap/CmdLine.h>
TCLAP::UnlabeledValueArg<std::string> dirArg("a", "dir", true, "models\\", "directory name");
TCLAP::UnlabeledValueArg<std::string> fileArg("b", "file", true, "filename", "file name");
TCLAP::ValueArg<std::string>          fieldArg("f", "field", "", false, "", "field path");
TCLAP::ValueArg<int>                  errorArg("e", "error", "compute error", false, 0, "int");
TCLAP::ValueArg<float>                thickArg("l", "", "base layer thickess", false, 0.6, "float");

TCLAP::SwitchArg screenshotArg("s", "screenshot", "generate screenshot", false);
TCLAP::SwitchArg deformedArg("d", "deformed", "", false);
TCLAP::SwitchArg reversededArg("r", "reversed", "", false);
TCLAP::ValueArg<float> yPlaneArg("y", "slice", "", false, 0.5, "float");
TCLAP::ValueArg<std::string>    posArg("p", "position", "screenshot position (F1, F2, ...), default: 1", false, "001", "int");

// --------------------------------------------------------------

#define SCREEN_W 1024
#define SCREEN_H 1024

AutoBindShader::shade     g_Sh;
GLuint g_ShAttribH;
GLuint g_ShAttribW;
GLuint g_ShAttribGradHdZ;

AutoBindShader::slicerror g_ShError;
GLuint g_ShErrorAttribH;
GLuint g_ShErrorAttribGradHdZ;

#include "TetMesh.h"

MeshFormat_stl                                 g_stlformat;
AutoPtr<TetMesh>                               g_Mesh;
//TriangleMesh_Ptr                               g_MeshStl;
//MeshRenderer<MeshFormat_stl::t_VertexFormat>*  g_Renderer;

float slice = 0.5f;

Array<float> h;

int displaymode = 1;
bool g_DrawTetsOrMesh = true;
bool g_DrawDeformed = false;
bool g_Screenshot = false;

Array<Tuple<double, 9>> g_grad_mats;

int   g_i_bottom_shape = -1;
float g_bottomz        = std::numeric_limits<float>::max();
float g_topz           = - std::numeric_limits<float>::max();
float g_extent_h       = 0.0f;

bool g_ComputeError = false;
enum e_ErrorType {e_Curved = 0, e_Adaptive = 1, e_Uniform = 2};
e_ErrorType g_ErrorToCompute = e_Curved;

// --------------------------------------------------------------

bool g_doScreenshot = false;
bool g_Reversed = false;

// --------------------------------------------------------------

map<v3u, int> tris_to_tets;

void buildTrisToTet()
{
  if (tris_to_tets.empty()) {
    ForIndex(t, g_Mesh->numTetrahedrons()) {
      if (!g_Mesh->isInside(t)) continue;
      ForIndex(off, 4) {
        v3u tri;
        ForIndex(i, 3) {
          tri[i] = g_Mesh->tetrahedronAt(t)[(off + i) % 4];
        }
        // sort tri
        if (tri[0] > tri[1]) std::swap(tri[0], tri[1]);
        if (tri[1] > tri[2]) std::swap(tri[1], tri[2]);
        if (tri[0] > tri[1]) std::swap(tri[0], tri[1]);
        tris_to_tets[tri] = t;
      }
    }
  }
}

// --------------------------------------------------------------

void mainKeyboard(unsigned char key)
{
  static float speed = 1.0f;
  static char last = ' ';

  if (key == 'q') {
    TrackballUI::exit();
  } else if (key == ' ') {
    displaymode = (displaymode + 1) % 4;
  } else if (key == 'm') {
    g_DrawTetsOrMesh = !g_DrawTetsOrMesh;
  } else if (key == 'd') {
    g_DrawDeformed = !g_DrawDeformed;
  } else if (key == 'r') {
    g_Reversed = !g_Reversed;
  } else if (key == 's') {
    g_Screenshot = true;
  } else if (key == '*') {
    g_ComputeError = !g_ComputeError;
  } else if (key == '/') {
    g_ErrorToCompute = (e_ErrorType)(((int)g_ErrorToCompute + 1) % 3);
  }

  if (key == '+') {
    slice = min(slice + 0.01, 1.00);
  }
  if (key == '-') {
    slice = max(slice - 0.01, 0.00);
  }

}

// --------------------------------------------------------------

void mainAnimate(double time, float elapsed)
{
}

// --------------------------------------------------------------

void mainOnReshape(uint x, uint y)
{

}

// --------------------------------------------------------------

void drawTriangle(v4f *pts, int a, int b, int c)
{
  glVertexAttrib1f(g_ShAttribH, pts[a][3]);
  glVertex3fv(&pts[a][0]);
  glVertexAttrib1f(g_ShAttribH, pts[b][3]);
  glVertex3fv(&pts[b][0]);
  glVertexAttrib1f(g_ShAttribH, pts[c][3]);
  glVertex3fv(&pts[c][0]);
}

// --------------------------------------------------------------

void render_tets()
{
  if (g_Reversed) {
    g_Sh.u_ClipDir.set(-1);
  } else {
    g_Sh.u_ClipDir.set(1);
  }

  float plane_y = slice * g_Mesh->getBBox().extent()[1] + g_Mesh->getBBox().minCorner()[1];

  glEnable(GL_CLIP_DISTANCE0);

  glDisable(GL_CULL_FACE);

  g_Sh.u_ClipY.set(plane_y + (g_Reversed ? 0.01f : -0.01f));
  glBegin(GL_TRIANGLES);
  ForIndex(t, g_Mesh->numTetrahedrons()) {
    ForIndex(off, 4) {
      ForIndex(i, 3) {
        int v = g_Mesh->tetrahedronAt(t)[(off + i) % 4];
        v3f pos = g_Mesh->vertexAt(v);
        glVertexAttrib1f(g_ShAttribH, h[v]);
        glVertexAttrib1f(g_ShAttribW, g_Mesh->isInside(t) ? 0.0f : 1.0f);
        glVertex3f(pos[0], pos[1], pos[2]);
      }
    }
  }
  glEnd();
  glDisable(GL_CULL_FACE);
  g_Sh.u_ClipY.set(plane_y + (g_Reversed ? -0.01f : 0.01f));
  g_Sh.u_Deform.set(g_DrawDeformed ? 1 : 0);
  glBegin(GL_TRIANGLES);
  ForIndex(t, g_Mesh->numTetrahedrons()) {

    AAB<3> tetbx;
    ForIndex(off, 4) {
      v3f pos = g_Mesh->vertexAt(g_Mesh->tetrahedronAt(t)[off]);
      tetbx.addPoint(pos);
    }

    // if outside and hits plane
    if (!g_Mesh->isInside(t)) {
      if (tetbx.minCorner()[1] <= plane_y && plane_y <= tetbx.maxCorner()[1]) {
        // intersect!
        int ni = 0;
        v4f intersects[16];
        ForIndex(o0, 4) {
          ForRange(o1, o0 + 1, 3) {
            v3f v0 = g_Mesh->vertexAt(g_Mesh->tetrahedronAt(t)[o0]);
            v3f v1 = g_Mesh->vertexAt(g_Mesh->tetrahedronAt(t)[o1]);
            v4f p0 = v4f(v0, h[g_Mesh->tetrahedronAt(t)[o0]]);
            v4f p1 = v4f(v1, h[g_Mesh->tetrahedronAt(t)[o1]]);
            if (p0[1] > p1[1]) {
              std::swap(p0, p1);
            }
            if (p0[1] < plane_y && plane_y <= p1[1]) {
              // intersection exists
              intersects[ni++] = (p0 + (p1 - p0) * ((plane_y - p0[1]) / (p1[1] - p0[1])));
              sl_assert(ni < 16);
            }
          }
        }

        glVertexAttrib1f(g_ShAttribW, 1.0f);
        if (ni == 3) {
          // draw!
          drawTriangle(intersects, 0, 1, 2);
        } else if (ni == 4) {
          // draw!
          drawTriangle(intersects, 0, 1, 2);
          drawTriangle(intersects, 0, 1, 3);
          drawTriangle(intersects, 0, 2, 3);
        }

      }
    }
  }
  glEnd();
  LIBSL_GL_CHECK_ERROR;
}

// --------------------------------------------------------------

void render_mesh()
{
  buildTrisToTet();

  // draw mesh
  glBegin(GL_TRIANGLES);
  glVertexAttrib1f(g_ShAttribW, 0.0f);
  ForIndex(i_srf, g_Mesh->numSurfaces()) {
    ForIndex(i, 3) {

      v3u tri = g_Mesh->surfaceTriangleAt(i_srf);
      v3u key = tri;
      if (key[0] > key[1]) std::swap(key[0], key[1]);
      if (key[1] > key[2]) std::swap(key[1], key[2]);
      if (key[0] > key[1]) std::swap(key[0], key[1]);

      int t = tris_to_tets[key];

      v4u tet = g_Mesh->tetrahedronAt(t);

      double ha = h[tet[0]];
      double hb = h[tet[1]];
      double hc = h[tet[2]];
      double hd = h[tet[3]];

      double tri_grad_h_dz =
        (ha - hd) * (float)g_grad_mats[t][6]
        + (hb - hd) * (float)g_grad_mats[t][7]
        + (hc - hd) * (float)g_grad_mats[t][8];

      v3f pos = g_Mesh->vertexAt(tri[i]);
      if (g_ErrorToCompute == e_Curved) {
        glVertexAttrib1f(g_ShAttribH, h[g_Mesh->surfaceTriangleAt(i_srf)[i]]);
        glVertexAttrib1f(g_ShAttribGradHdZ, 1.0f / (float)tri_grad_h_dz);
      } else {
        glVertexAttrib1f(g_ShAttribGradHdZ, 1.0f);
        glVertexAttrib1f(g_ShAttribH, pos[2]);
      }

      if (g_DrawDeformed) {
        glVertex3f(pos[0], pos[1], pos[2]);
      } else {
        glVertex3f(pos[0], pos[1], h[g_Mesh->surfaceTriangleAt(i_srf)[i]]);
      }
    }
  }
  glEnd();
}

// --------------------------------------------------------------

void makeAdaptiveSequence(int num_layers, Array<v3d>& sequence, GLTexBuffer& tex_sequence)
{
  // static int num_layers = 0;
  // num_layers = (num_layers+1)%155;
  if (g_ErrorToCompute == e_Adaptive) {
    string filename = dirArg.getValue() + fileArg.getValue() + sprint("/%.1f/sequence_%04d.array_v3d",thickness,num_layers);
    if (!LibSL::System::File::exists(filename.c_str())) {
      cerr << "cannot find sequence " << filename << endl;
      exit(-1);
    }
    loadArray(sequence, filename.c_str());
  }
  if (!sequence.empty()) {
    tex_sequence.init(2 * sequence.size() * sizeof(float), GL_R32F);
    Array<v2f> rawseq(sequence.size());
    ForIndex(i, sequence.size()) {
      rawseq[i][0] = (float)sequence[i][0];
      rawseq[i][1] = (float)sequence[i][1];
    }
    tex_sequence.writeTo(rawseq);
  }
}

// --------------------------------------------------------------

// This below is highly innefficient, but that is not the goal is it?

double render_error()
{
 // #define DEBUG_ERROR

  GLProtectViewport vp;

  RenderTarget2DRGBA rt(2048,2048);

  int num_layers = ceil(g_extent_h / thickness);

  Array<v3d>  sequence;
  GLTexBuffer tex_sequence;
  if (g_ErrorToCompute == e_Adaptive) {
    makeAdaptiveSequence(num_layers, sequence, tex_sequence);
  }

#ifndef DEBUG_ERROR
  rt.bind();
  glViewport(0, 0, rt.w(), rt.h());
#endif

  glClearColor(0, 0, 0, 0);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE);

  GLTexBuffer  errors;
  errors.init(rt.w()*rt.h()*sizeof(uint), GL_R32UI);
  Array<uint> rawdata(rt.w()*rt.h());
  rawdata.fill(0);
  errors.writeTo(rawdata);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  AAB<3> bx = g_Mesh->getBBox();

  g_ShError.begin();

  // SL: NOTE: view is distorted but this is taken into account below when computing the final error
  // float max_ex = max(bx.extent()[0], bx.extent()[1]);
  float l = bx.minCorner()[0];//bx.center()[0] - max_ex / 2.0f;
  float r = bx.maxCorner()[0];//bx.center()[0] + max_ex / 2.0f;
  float t = bx.minCorner()[1];//bx.center()[1] - max_ex / 2.0f;
  float b = bx.maxCorner()[1];//bx.center()[1] + max_ex / 2.0f;

  g_ShError.u_matrix.set(
    orthoMatrixGL(
      l,r,t,b,-bx.minCorner()[2], -bx.maxCorner()[2]
    )
  );

  //g_ShError.u_matrix.set(perspectiveMatrixGL(float(M_PI / 8.0), 1.0f, 1.0f, 1000.0f)    * TrackballUI::matrix());

  if (g_ErrorToCompute == e_Curved) {
    g_ShError.u_v_h_min.set(h[g_i_bottom_shape]);
  } else {
    g_ShError.u_v_h_min.set(g_bottomz);
  }

  float thu = (g_topz - g_bottomz) / (float)num_layers;
  if (g_ErrorToCompute == e_Uniform) {
    // compute thickness giving the correct number of layers
    g_ShError.u_thickness.set( thu );
  } else {
    g_ShError.u_thickness.set( (float)thickness);
  }
  g_ShError.u_base_thickness.set((float)thickness);

  glBindImageTexture(0, errors.glTexId(), 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);
  g_ShError.u_Errors.set(0);

  if (!sequence.empty()) {
    glBindImageTexture(1, tex_sequence.glTexId(), 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);
    g_ShError.u_Sequence.set(1);
    g_ShError.u_SequenceLength.set((int)sequence.size());
  } else {
    g_ShError.u_SequenceLength.set(0);
  }

  // draw mesh
  buildTrisToTet();

  glBegin(GL_TRIANGLES);
  ForIndex(i_srf, g_Mesh->numSurfaces()) {

    v3u tri = g_Mesh->surfaceTriangleAt(i_srf);

    v3u key = tri;
    if (key[0] > key[1]) std::swap(key[0], key[1]);
    if (key[1] > key[2]) std::swap(key[1], key[2]);
    if (key[0] > key[1]) std::swap(key[0], key[1]);

    int t = tris_to_tets[key];

    v4u tet = g_Mesh->tetrahedronAt(t);

    double ha = h[tet[0]];
    double hb = h[tet[1]];
    double hc = h[tet[2]];
    double hd = h[tet[3]];

    double tri_grad_h_dz =
        (ha - hd) * (float)g_grad_mats[t][6]
      + (hb - hd) * (float)g_grad_mats[t][7]
      + (hc - hd) * (float)g_grad_mats[t][8];

    ForIndex(i, 3) {
      v3f pos = g_Mesh->vertexAt(tri[i]);

      if (g_ErrorToCompute == e_Curved) {
        glVertexAttrib1f(g_ShErrorAttribGradHdZ, 1.0f / (float)tri_grad_h_dz);
        glVertexAttrib1f(g_ShErrorAttribH, h[g_Mesh->surfaceTriangleAt(i_srf)[i]]);
      } else {
        glVertexAttrib1f(g_ShErrorAttribGradHdZ, 1.0f);
        glVertexAttrib1f(g_ShErrorAttribH, pos[2]);
      }
      glVertex3f(pos[0], pos[1], pos[2]);
    }
  }
  glEnd();

  g_ShError.end();

#ifndef DEBUG_ERROR
  rt.unbind();
#endif

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  glEnable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

#ifndef DEBUG_ERROR
  errors.readBack(rawdata);

  errors.terminate();

  const double fxp = 4096.0;
  double sum_of_all = 0.0;

  ForArray(rawdata, i) {
    sum_of_all += thickness * 0.5 * ((double)rawdata[i]) / fxp;
  }

  double volume_error = sum_of_all * (bx.extent()[0] * bx.extent()[1]) / (rt.w() * rt.h());

  if (g_ErrorToCompute == e_Uniform)  cerr << "uniform  ";
  if (g_ErrorToCompute == e_Curved)   cerr << "curved   ";
  if (g_ErrorToCompute == e_Adaptive) cerr << "adaptive ";
  cerr << " error : " << volume_error << " mm^3";
  if (g_ErrorToCompute == e_Uniform)  cerr << " ( num layers : " << ceil((g_topz - g_bottomz) / thu) << ", thickness : " << thu << ") ";
  if (g_ErrorToCompute == e_Adaptive) cerr << " ( num layers : " << num_layers << " )";
  cerr << endl;

  // store error in file
  if (errorArg.isSet()) {
    std::string freport = dirArg.getValue() + "..\\" + fileArg.getValue() + ".txt";
    ofstream f(freport, ios::app);
    if (g_ErrorToCompute == e_Curved) {
      // external script is expected to call on curve first, so this make a header
      f << endl;
      f << fieldArg.getValue() << endl;
    }
    if (g_ErrorToCompute == e_Uniform)  f << "uniform  ";
    if (g_ErrorToCompute == e_Curved)   f << "curved   ";
    if (g_ErrorToCompute == e_Adaptive) f << "adaptive ";
    f << " error : " << volume_error << " mm^3";
    if (g_ErrorToCompute == e_Uniform)  f << " ( num layers : " << ceil((g_topz - g_bottomz) / thu) << ", thickness : " << thu << ") ";
    if (g_ErrorToCompute == e_Adaptive) f << " ( num layers : " << num_layers << " )";
    f << " [" << thickness << "] ";
    f << endl;
    f.close();
    exit(0);
  }

  return volume_error;
#else
  return 0.0f;
#endif
}

// --------------------------------------------------------------

void render_bbox()
{
  glDisable(GL_CLIP_DISTANCE0);

  /// draw bbox
  v3f min = g_Mesh->getBBox().minCorner();
  v3f max = g_Mesh->getBBox().maxCorner();

  glLineWidth(2.5);
  glColor3f(1.0, 1.0, 1.0);

  glBegin(GL_LINES);

  glVertex3f(min[0], min[1], min[2]); glVertex3f(max[0], min[1], min[2]);
  glVertex3f(max[0], min[1], min[2]); glVertex3f(max[0], max[1], min[2]);
  glVertex3f(max[0], max[1], min[2]); glVertex3f(min[0], max[1], min[2]);
  glVertex3f(min[0], max[1], min[2]); glVertex3f(min[0], min[1], min[2]);

  glVertex3f(min[0], min[1], max[2]); glVertex3f(max[0], min[1], max[2]);
  glVertex3f(max[0], min[1], max[2]); glVertex3f(max[0], max[1], max[2]);
  glVertex3f(max[0], max[1], max[2]); glVertex3f(min[0], max[1], max[2]);
  glVertex3f(min[0], max[1], max[2]); glVertex3f(min[0], min[1], max[2]);

  glVertex3f(min[0], min[1], min[2]); glVertex3f(min[0], min[1], max[2]);
  glVertex3f(max[0], min[1], min[2]); glVertex3f(max[0], min[1], max[2]);
  glVertex3f(max[0], max[1], min[2]); glVertex3f(max[0], max[1], max[2]);
  glVertex3f(min[0], max[1], min[2]); glVertex3f(min[0], max[1], max[2]);

  glEnd();
}

// --------------------------------------------------------------

void take_screenshot()
{
  ImageRGBA_Ptr img = ImageRGBA_Ptr(new ImageRGBA(SCREEN_W, SCREEN_H));
  glReadPixels(0, 0, SCREEN_W, SCREEN_H, GL_RGBA, GL_UNSIGNED_BYTE, img->pixels().raw());

  ImageRGB_Ptr rgb = ImageRGB_Ptr(img->cast<ImageRGB>());
  rgb->flipH();

  static int cnt = 0;
  string file;
  while (cnt < 4096) {
    file = sprint("shot%04d.png", cnt++);
    if (!LibSL::System::File::exists(file.c_str())) {
      break;
    }
  }
  saveImage(file.c_str(), rgb);
  cerr << "Screenshot saved as: " << file << endl;

  if (g_doScreenshot) {
    exit(0);
  }
}

// --------------------------------------------------------------

void mainRender()
{
  /// render on screen

  LIBSL_GL_CHECK_ERROR;

  GPUHelpers::clearScreen(LIBSL_COLOR_BUFFER | LIBSL_DEPTH_BUFFER,
    1, 1, 1);

  if (g_ComputeError) {
    render_error();
    return;
  }

  LIBSL_GL_CHECK_ERROR;

  m4x4f mat =
    perspectiveMatrixGL(float(M_PI / 8.0), 1.0f, 1.0f, 1000.0f)
    * TrackballUI::matrix()
    * translationMatrix(-g_Mesh->getBBox().center())
    ;
  v4f eye_at = mat.inverse() * v4f(0, 0, 1, 0);
  float dist = length(eye_at - mat.inverse() * v4f(0, 0, -1, 0));

  LIBSL_GL_CHECK_ERROR;

  g_Sh.begin();
  g_Sh.u_model.set(TrackballUI::matrix());
  g_Sh.u_matrix.set(mat);
  float plane_y = slice * g_Mesh->getBBox().extent()[1] + g_Mesh->getBBox().minCorner()[1];

  g_Sh.u_drawing_mesh.set(g_DrawTetsOrMesh ? 0 : 1);

  LIBSL_GL_CHECK_ERROR;

  /// draw tets

  glDisable(GL_CLIP_DISTANCE0);

  if (g_DrawTetsOrMesh) {

    render_tets();
    LIBSL_GL_CHECK_ERROR;

  } else {

    if (g_ErrorToCompute == e_Curved) {
      g_Sh.u_v_h_min.set(h[g_i_bottom_shape]);
    } else {
      g_Sh.u_v_h_min.set(g_bottomz);
    }
    int num_layers = ceil(g_extent_h / thickness);
    if (g_ErrorToCompute == e_Uniform) {
      // compute thickness giving the correct number of layers
      float thu = (g_topz - g_bottomz) / (float)num_layers;
      g_Sh.u_thickness.set(thu);
    } else {
      g_Sh.u_thickness.set((float)thickness);
    }
    g_Sh.u_base_thickness.set((float)thickness);

    Array<v3d>  sequence;
    GLTexBuffer tex_sequence;
    if (g_ErrorToCompute == e_Adaptive) {
      makeAdaptiveSequence(num_layers, sequence, tex_sequence);
    }
    if (!sequence.empty()) {
      glBindImageTexture(1, tex_sequence.glTexId(), 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);
      g_Sh.u_Sequence.set(1);
      g_Sh.u_SequenceLength.set((int)sequence.size());
    } else {
      g_Sh.u_SequenceLength.set(0);
    }

    render_mesh();
    LIBSL_GL_CHECK_ERROR;

  }

  if (g_DrawTetsOrMesh && !g_DrawDeformed) {

    render_bbox();
    LIBSL_GL_CHECK_ERROR;
  }

  g_Sh.end();

  if (g_Screenshot) {

    take_screenshot();

    g_Screenshot = false;
  }

  LIBSL_GL_CHECK_ERROR;

}

/* -------------------------------------------------------- */

int main(int argc, char **argv)
{
  std::string additionnal = "";

  TCLAP::CmdLine cmd("", ' ', "1.0");

  cmd.add(dirArg);
  cmd.add(fileArg);
  cmd.add(fieldArg);
  cmd.add(errorArg);
  cmd.add(thickArg);

  cmd.add(screenshotArg);
  cmd.add(deformedArg);
  cmd.add(reversededArg);
  cmd.add(yPlaneArg);
  cmd.add(posArg);

  cmd.parse(argc, argv);

  g_Screenshot   = g_doScreenshot = screenshotArg.isSet();
  g_DrawDeformed = deformedArg.isSet();
  g_Reversed     = reversededArg.isSet();
  slice = yPlaneArg.isSet() ? yPlaneArg.getValue() : 0.5;
  thickness = thickArg.isSet() ? thickArg.getValue() : 0.6;

  try {

    /// init TrackballUI UI (glut clone for both GL and D3D, with a trackball)
    cerr << "Init TrackballUI   ";
    TrackballUI::onRender = mainRender;
    TrackballUI::onKeyPressed = mainKeyboard;
    TrackballUI::onAnimate = mainAnimate;
    TrackballUI::onReshape = mainOnReshape;
    TrackballUI::init(SCREEN_W, SCREEN_H);

    if (posArg.isSet()) {
      TrackballUI::trackballLoad( posArg.getValue().c_str() );
    }

    cerr << "[OK]" << endl;

    /// help
    printf("[q]     - quit\n");

    glEnable(GL_DEPTH_TEST);


    g_ShError.init();
    g_ShErrorAttribH = glGetAttribLocationARB(g_ShError.shader().handle(), "a_h");
    g_ShErrorAttribGradHdZ = glGetAttribLocationARB(g_ShError.shader().handle(), "a_grad_h_z");

    g_Sh.init();
    g_ShAttribH = glGetAttribLocationARB(g_Sh.shader().handle(), "a_h");
    g_ShAttribW = glGetAttribLocationARB(g_Sh.shader().handle(), "a_widding");
    g_ShAttribGradHdZ = glGetAttribLocationARB(g_ShError.shader().handle(), "a_grad_h_z");

    LIBSL_GL_CHECK_ERROR;

    // open mesh
    cerr << "Loading mesh      ";

    string filename = dirArg.getValue() + fileArg.getValue() + ".msh";
    cerr << filename << " ";

    g_Mesh = AutoPtr<TetMesh>(new TetMesh());
    g_Mesh = AutoPtr<TetMesh>(g_Mesh->load(filename.c_str()));

    cerr << "[OK]" << endl;
    cerr << sprint("  mesh contains %d vertices, %d triangles, %d tets\n", g_Mesh->numVertices(), g_Mesh->numTriangles(), g_Mesh->numTetrahedrons());

    cerr << "bbox center: " << g_Mesh->getBBox().center() << endl;
    cerr << "bbox min corner: " << g_Mesh->getBBox().minCorner() << endl;
    cerr << "bbox max corner: " << g_Mesh->getBBox().maxCorner() << endl;

    cerr << "Creating renderer " << endl;

    string displacement = fieldArg.isSet()
      ? fieldArg.getValue()
      : dirArg.getValue() + fileArg.getValue() +  "/displacements" ;
    loadArray(h, displacement.c_str() );

    string tetmats = dirArg.getValue() + fileArg.getValue() + "/tetmats";

    loadArray(g_grad_mats, tetmats.c_str());

    // i_bottom_shape is the variable at the bottom of the shape
    // this is needed for slice alignment
    float hmin_shape = std::numeric_limits<float>::max();
    float hmax_shape = -std::numeric_limits<float>::max();
    ForIndex(t, g_Mesh->numSurfaces()) {
      v3u tri = g_Mesh->surfaceTriangleAt(t);
      ForIndex(i, 3) {
        v3f p = g_Mesh->vertexAt(tri[i]);
        if (p[2] < g_bottomz) {
          g_bottomz = p[2];
          g_i_bottom_shape = tri[i];
        }
        g_topz = max(p[2], g_topz);
        hmin_shape = min(hmin_shape, h[tri[i]]);
        hmax_shape = max(hmax_shape, h[tri[i]]);
      }
    }
    g_extent_h = (hmax_shape - hmin_shape);
    cerr << "model z extent in deformed space: " << g_extent_h << endl;

    g_Sh.begin();
    g_Sh.u_v_h_min.set(h[g_i_bottom_shape]);
    g_Sh.end();

    cerr << "[OK]" << endl;

    // setup view

    float h_avg = 0.0f;
    ForIndex(i, h.size()) {
      h_avg += h[i];
    }
    h_avg /= (float)h.size();

    // init trackball viewpoint
    TrackballUI::trackball().setCenter(
      //v3f(v2f(g_Mesh->getBBox().center()), h_avg)
      v3f(g_Mesh->getBBox().center())
    );

    if (errorArg.isSet()) {
      g_ComputeError = true;
      g_ErrorToCompute = (e_ErrorType)errorArg.getValue();
    }

    /// main loop
    TrackballUI::loop();

    /// clean exit
    g_Sh.terminate();
    g_ShError.terminate();

    // shutdown UI
    TrackballUI::shutdown();

  } catch (Fatal& e) {
    cerr << e.message() << endl;
    return (-1);
  }

  return (0);
}

/* -------------------------------------------------------- */
