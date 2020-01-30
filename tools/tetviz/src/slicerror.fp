#version 430 core
#extension GL_ARB_shader_image_load_store : enable
#extension GL_ARB_gpu_shader5 : enable

uniform uint  u_ScreenW;
uniform float u_v_h_min;
uniform int   u_SequenceLength;

layout(r32ui) coherent uniform uimageBuffer u_Errors;
layout(r32f)  coherent uniform imageBuffer  u_Sequence;

in float   v_h;
in float   v_grad_h_z;

uniform float u_base_thickness;
uniform float u_thickness;

vec2 findSlice(float z)
{
  int l = 0;
  int r = u_SequenceLength;
  for (int i = 0 ; i < 64 ; i++) {
    // read sequence
	int m = min((l+r)/2,u_SequenceLength-1);
	vec2 slice_nfo;
	slice_nfo.x = imageLoad(u_Sequence,m*2+0).x;
	slice_nfo.y = imageLoad(u_Sequence,m*2+1).x;
	if (z >= slice_nfo.x && z < slice_nfo.x + slice_nfo.y) {
	  return slice_nfo;
	}
	if (z < slice_nfo.x) {
	  r = m;
	} else {
	  l = m;
	}
  }
  int last = u_SequenceLength-1;
  return vec2(imageLoad(u_Sequence,last*2+0).x,imageLoad(u_Sequence,last*2+1).x);
}

void main()
{
  float normalized_distance =  0.0;

  if (u_SequenceLength > 0) {
    
	vec2 slice = findSlice( v_h );
    normalized_distance = (slice.y / u_base_thickness) * (1.0 - abs(0.5 - ( (v_h-slice.x) / slice.y) ) * 2.0);

  } else {

    normalized_distance = v_grad_h_z * (u_thickness / u_base_thickness) * (1.0 - abs(0.5 - fract( (v_h-u_v_h_min) / u_thickness) ) * 2.0);

  }
  
  uint  ptr = uint(gl_FragCoord.x) + uint(gl_FragCoord.y) * u_ScreenW;
  uint  fxp_distance = uint(normalized_distance * 4096.0);
  imageAtomicAdd( u_Errors, int(ptr), fxp_distance );

  gl_FragColor = vec4( normalized_distance );
}
