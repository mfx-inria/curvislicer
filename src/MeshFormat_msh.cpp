//---------------------------------------------------------------------------
#include "LibSL.precompiled.h"
//---------------------------------------------------------------------------

#include "MeshFormat_msh.h"
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

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Mesh

//---------------------------------------------------------------------------

/// Declaring a global will automatically register the plugin
namespace {
  NAMESPACE::MeshFormat_msh s_Msh;  /// FIXME: this mechanism does not work with VC++
}                                    ///        see also MeshFormatManager constructor

//---------------------------------------------------------------------------

NAMESPACE::MeshFormat_msh::MeshFormat_msh()
{
  try {
    // register plugin
    MESH_FORMAT_MANAGER.registerPlugin(this);
  } catch (LibSL::Errors::Fatal& e) {
    std::cerr << e.message() << std::endl;
  }
}

//---------------------------------------------------------------------------

using namespace std;

//---------------------------------------------------------------------------

void parseWhiteSpace(ifstream& stream) {
  char next = stream.peek();
  while (next == '\n' || next == ' ' || next == '\t' || next == '\r') {
    stream.get();
    next = stream.peek();
  }
}

void parseMeshFormat(ifstream& stream, MeshFormat_msh::msh_format_info& format) {
  stream >> format.version_number;
  parseWhiteSpace(stream);
  stream >> format.binary;
  parseWhiteSpace(stream);
  stream >> format.data_size;
  parseWhiteSpace(stream);

  if (format.binary) { // binary
    int check = 1;
    stream.read(reinterpret_cast<char*>(&check), sizeof(int));
    sl_assert(check == 1); // make sure the file is compatible (little or big endian)
  }
}

void parsePhysicalNames(ifstream& stream, MeshFormat_msh::msh_format_info& format)
{
}

void parseNodes(ifstream& stream, MeshFormat_msh::msh_format_info& format)
{
  uint nb_nodes;
  stream >> nb_nodes;
  parseWhiteSpace(stream);

  int num_i;
  double xyz[3];
  
  if (format.binary) {
    ForIndex(i, nb_nodes) {
	    stream.read(reinterpret_cast<char*>(&num_i), sizeof(int));
	    stream.read(reinterpret_cast<char*>(&xyz),   sizeof(xyz));
	    format.nodes.emplace_back(v3f(xyz[0], xyz[1], xyz[2]));
    }
  } else {
	  ForIndex (i, nb_nodes) {
      char trash;
      stream >> num_i;
	    stream >> xyz[0] >> trash >> xyz[1] >> trash >> xyz[2];
	    parseWhiteSpace(stream);
	  }
  }
}

void parseElements(ifstream& stream, MeshFormat_msh::msh_format_info& format)
{
  uint nb_elements;
  stream >> nb_elements;
  parseWhiteSpace(stream);

  int header[3];
  uint nb_elems = 0;
  while (nb_elems < nb_elements)
  {
    
	  if (format.binary) {
	    stream.read(reinterpret_cast<char*>(&header), sizeof(header));
    } else {
	    stream >> header[0];
	    parseWhiteSpace(stream);
	  }
	
	  int elem_size = 0;
	
    switch (header[0]) {
    case 1: // 2-node line
      break;
    case 2: // 3-node triangle
	    elem_size = 3;
      break;
    case 3: // 4-node quadrangle
      elem_size = 4;
	  break;
    case 4: // 4-node tetrahedron 
      elem_size = 4;
      break;
    case 5: // 8-node hexahedron
      break;
    case 6: 
      break;
    }
	
	for (int elem = 0; elem < header[1]; elem++) {
      for (int tag = 0; tag < header[2]; tag++) {
        int data;
        stream.read(reinterpret_cast<char*>(&data), sizeof(data));
      }
      int num_i;
      stream.read(reinterpret_cast<char*>(&num_i), sizeof(num_i));

      int data[4];
      stream.read(reinterpret_cast<char*>(&data), sizeof(data));
      vector<int> elems;
	  
	    ForIndex (i, elem_size) {
	      elems.emplace_back(data[i]);
	    }
	    format.elements.emplace_back(elems);
    }
    nb_elems += header[1];
  }
}

void parseElementData(ifstream& stream, MeshFormat_msh::msh_format_info& format) {
  MeshFormat_msh::msh_entity_data element;

  {
    uint nb_string_tags;
    stream >> nb_string_tags;
    std::string tag;
    ForIndex(i, nb_string_tags) {
      string str;
      do {
        stream >> str;
        tag += (tag.empty() ? "" : " ") + str;
      } while (str.at(str.size() - 1) != '\"');
      ;
      element.string_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }
  {
    uint nb_real_tags;
    stream >> nb_real_tags;
    float tag;
    ForIndex(i, nb_real_tags) {
      stream >> tag;
      element.real_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }
  {
    uint nb_int_tags;
    stream >> nb_int_tags;
    int tag;
    ForIndex(i, nb_int_tags) {
      stream >> tag;
      element.integer_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }

  size_t nb_values = element.integer_tags[1];
  size_t nb_elems  = element.integer_tags[2];

  size_t num_bytes = (nb_values * format.data_size + 4) * nb_elems;
  char* data = new char[num_bytes];
  stream.read(data, num_bytes);


  ForIndex(i, nb_elems) {
    size_t base_idx = i * (4 + nb_values * format.data_size);
    int elem = (*reinterpret_cast<int*>(&data[base_idx])) - 1;
    base_idx += 4;

    vector<float> values(nb_values);
    ForIndex(j, nb_values) {
      values[j] = (*reinterpret_cast<double*>(&data[base_idx + j * format.data_size]));
    }
    element.values.emplace(elem, values);
  }
  format.element_data.emplace_back(element);
}




void parseNodeData(ifstream& stream, MeshFormat_msh::msh_format_info& format) {
  MeshFormat_msh::msh_entity_data node_data;

  {
    uint nb_string_tags;
    stream >> nb_string_tags;
    std::string tag;
    ForIndex(i, nb_string_tags) {
      string str;
      do {
        stream >> str;
        tag += (tag.empty() ? "" : " ") + str;
      } while (str.at(str.size() - 1) != '\"');
     ;
      node_data.string_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }
  {
    uint nb_real_tags;
    stream >> nb_real_tags;
    float tag;
    ForIndex(i, nb_real_tags) {
      stream >> tag;
      node_data.real_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }
  {
    uint nb_int_tags;
    stream >> nb_int_tags;
    int tag;
    ForIndex(i, nb_int_tags) {
      stream >> tag;
      node_data.integer_tags.emplace_back(tag);
      parseWhiteSpace(stream);
    }
  }

  size_t nb_values = node_data.integer_tags[1];
  size_t nb_nodes  = node_data.integer_tags[2];

  size_t num_bytes = (nb_values * format.data_size + 4) * nb_nodes;
  char* data = new char[num_bytes];
  stream.read(data, num_bytes);

  ForIndex (i, nb_nodes) {
    int node = *reinterpret_cast<int*>(&data[i*(4 + nb_values * format.data_size)]);
    node -= 1;
    size_t base_idx = i * (4 + nb_values * format.data_size) + 4;

    vector<float> values(nb_values);
    ForIndex (j, nb_values) {
      values[j] = (*reinterpret_cast<float*>(&data[base_idx + j * format.data_size]));
    }

    node_data.values.emplace(node, values);
  }

  format.node_data.emplace_back(node_data);
}
//---------------------------------------------------------------------------

void NAMESPACE::MeshFormat_msh::parse(const char *fname)
{
  LIBSL_BEGIN;
  // open file
  FILE *f = NULL;
  fopen_s(&f, fname, "rb");
  if (f == NULL) {
    throw Fatal("[MeshFormat_msh::load] - file '%s' not found", fname);
  }

  ifstream stream(fname, ios::binary);

  std::string blockname, blockend;
  while (!stream.eof()) {
    stream >> blockname;

    if (blockname == "$MeshFormat") {
      parseMeshFormat(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndMeshFormat");
    } else if (blockname == "$PhysicalNames") {
      parsePhysicalNames(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndPhysicalNames");
    } else if (blockname == "$Nodes") {
      parseNodes(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndNodes");
    } else if (blockname == "$Elements") {
      parseElements(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndElements");
    } else if (blockname == "$ElementData") {
      parseElementData(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndElementData");
    } else if (blockname == "$NodeData") {
      parseNodeData(stream, m_format);
      stream >> blockend;
      sl_assert(blockend == "$EndNodeData");
    } else {
      cerr << "Unknown block : " << blockname << endl;
      // TODO: READ UNTIL END OF BLOCK
      stream >> blockend;
      sl_assert(blockend == "$End" + blockname.substr(1, blockname.size() - 1));
    }
    parseWhiteSpace(stream);
  }
  fclose(f);
  LIBSL_END;
}


NAMESPACE::TriangleMesh* NAMESPACE::MeshFormat_msh::load(const char *fname) const
{
  MeshFormat_msh msh;
  msh.parse(fname);

  msh_format_info format = msh.m_format;

  set<v3u> tris;
  map<v3u, uint> srfTris;

  // faces
  uint fi = 0;
  ForArray(format.elements, n_tet) { // TODO: manage non tet
    ForIndex(off, 4) {
      uint A = format.elements.at(n_tet)[(off + 0) % 4] - 1;
      uint B = format.elements.at(n_tet)[(off + 1) % 4] - 1;
      uint C = format.elements.at(n_tet)[(off + 2) % 4] - 1;

      if (A > B) std::swap(A, B);
      if (B > C) std::swap(B, C);
      if (A > B) std::swap(A, B);

      // extract the surface
      v3u t(A, B, C);

      auto it_ = tris.find(t);
      if (it_ == tris.end()) {
        tris.insert(t);

        fi++;
      }

      auto it = srfTris.find(t);
      if (it == srfTris.end()) {
        srfTris.insert_or_assign(t, fi);
      } else {
        srfTris.erase(it);
      }

    }
  }

  TriangleMesh_generic<MeshFormat_msh::t_VertexData> *mesh
    = new TriangleMesh_generic<MeshFormat_msh::t_VertexData>(format.nodes.size(), tris.size(), 1 , AutoPtr<MVF>(MVF::make<MeshFormat_msh::t_VertexFormat>()));

  //vertices
  ForArray (format.nodes, n_ver) {
    v3f pos = format.nodes.at(n_ver);
    mesh->vertexAt(n_ver).pos = pos;
  }

  // tris
  fi = 0;
  for (auto it : tris) {
    mesh->triangleAt(fi) = it;
    fi++;
  }

  //surface tris
  uint si = 0;
  TriangleMesh::t_SurfaceNfo srf;
  srf.triangleIds.allocate(srfTris.size());
  for(auto it : srfTris) {
    srf.triangleIds[si] = it.second;
    si++;
  }
  mesh->surfaceAt(0) = srf;

	// done
	return mesh;
}

//---------------------------------------------------------------------------

void NAMESPACE::MeshFormat_msh::save(const char *fname,const NAMESPACE::TriangleMesh *mesh) const
{
}

//---------------------------------------------------------------------------