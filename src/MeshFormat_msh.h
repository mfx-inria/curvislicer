// ------------------------------------------------------
// LibSL::Mesh::MeshFormat_msh
// ------------------------------------------------------
//
// Load/save LibSL 'msh' file format
//
// ------------------------------------------------------
// Jimmy ETIENNE - 2018-10-11
// ------------------------------------------------------

#pragma once

#include <LibSL/Math/Vertex.h>
#include <LibSL/Mesh/Mesh.h>
#include <LibSL/Math/Tuple.h>

namespace LibSL {
  namespace Mesh {

    class MeshFormat_msh : public TriangleMeshFormat_plugin
    {

    public:
      typedef struct
      {
        LibSL::Math::v3f pos;
        LibSL::Math::v3f nrm;
      } t_VertexData;

      typedef MVF2(mvf_position_3f, mvf_normal_3f) t_VertexFormat;

      typedef struct
      {
        std::vector<std::string> string_tags;  // [name of the post-processing view, name of the interpolation scheme, ???]
        std::vector<float>       real_tags;    // [time value, ???]
        std::vector<int>         integer_tags; // [time step index, number of field components (1, 3 or 9), number of entities, partition index, ???] 

        std::map<int, std::vector<float>> values; // entity number : [value1, value2, ...]

      } msh_entity_data;


      typedef struct
      {
        std::vector<std::string> string_tags;  // [name of the post-processing view, name of the interpolation scheme, ???]
        std::vector<float>       real_tags;    // [time value, ???]
        std::vector<int>         integer_tags; // [time step index, number of field components (1, 3 or 9), number of elements, partition index, ???] 

        std::map<int, std::pair<int, std::vector<float>>> values; // element number : <number of nodes per element ,[value1, value2, ...]>

      } msh_element_node_data;


      typedef struct
      {
        float  version_number;
        bool   binary; // 0 = ASCII, 1 = binary
        
        uint   data_size = 0;
        uint   number_of_nodes_per_element;
        uint   element_type;

        std::vector<LibSL::Math::v3f> nodes;
        std::vector<std::vector<int>> elements;

        std::vector<msh_entity_data>       node_data;
        std::vector<msh_entity_data>       element_data;
        std::vector<msh_element_node_data> element_node_data;

      } msh_format_info;

    public:
      MeshFormat_msh();
      void           parse(const char *);
      void           save(const char *,const TriangleMesh *) const;
      TriangleMesh  *load(const char *)                      const;
      const char    *signature()                             const {return "msh";}
    
    public:
      mutable std::vector<uint>                        m_meshes;
      mutable std::vector<LibSL::Math::Tuple<uint, 3>> m_surfaces;
      msh_format_info                                  m_format;
    };

  } //namespace LibSL::Mesh
} //namespace LibSL

// ------------------------------------------------------
