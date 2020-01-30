#version 430 core
#extension GL_ARB_shader_image_load_store : enable
#extension GL_ARB_gpu_shader5 : enable

in float   v_h;
in float   v_grad_h_z;

in float   v_z;

in vec3    v_pos;
in float   v_widding;

uniform int   u_drawing_mesh;
uniform float u_v_h_min;

uniform float u_base_thickness;
uniform float u_thickness;

uniform int   u_SequenceLength;
layout(r32f)  coherent uniform imageBuffer  u_Sequence;

vec2 findSlice(float z)
{
  int l = 0;
  int r = u_SequenceLength;
  for (int i = 0; i < 64; i++) {
    // read sequence
    int m = min((l + r) / 2, u_SequenceLength - 1);
    vec2 slice_nfo;
    slice_nfo.x = imageLoad(u_Sequence, m * 2 + 0).x;
    slice_nfo.y = imageLoad(u_Sequence, m * 2 + 1).x;
    if (z >= slice_nfo.x && z < slice_nfo.x + slice_nfo.y) {
      return slice_nfo;
    }
    if (z < slice_nfo.x) {
      r = m;
    } else {
      l = m;
    }
  }
  int last = u_SequenceLength - 1;
  return vec2(imageLoad(u_Sequence, last * 2 + 0).x, imageLoad(u_Sequence, last * 2 + 1).x);
}

void main()
{
  float t = fract(v_h / 5);

  vec3 dpdx = dFdx(v_pos);
  vec3 dpdy = dFdy(v_pos);
  vec3 nrm  = normalize( cross(dpdx,dpdy));

  if (u_drawing_mesh == 1) {

    if (u_SequenceLength > 0) {

      vec2 slice = findSlice(v_h);
      gl_FragColor = (slice.y / u_base_thickness) * vec4(1.0 - abs(0.5 - ((v_h - slice.x) / slice.y)) * 2.0);

    } else {

      gl_FragColor = v_grad_h_z * (u_thickness / u_base_thickness) * vec4(1.0 - abs(0.5 - fract((v_h - u_v_h_min) / u_thickness)) * 2.0);

    }

  } else {

  if (v_widding > 0.5) {
    if (t < 0.5) {
      gl_FragColor = vec4(0.8, 0.4, 0.6, 0)*(0.5+0.5*nrm.z);
    } else {
      gl_FragColor = vec4(0, 1, 1, 0)*(0.5 + 0.5*nrm.z);
    }
  } else {
    if (t < 0.5) {
      gl_FragColor = vec4(0.4, 0.8, 0.6, 0)*nrm.z;
    } else {
      gl_FragColor = vec4(0.2, 1, 0.6, 0)*nrm.z;
    }
  }

  }
}
