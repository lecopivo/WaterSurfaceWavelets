const int NUM = DIR_NUM/4;
const int NUM_INTEGRATION_NODES = 4*DIR_NUM;

// Define Ampl for vertex shader
#if VERTEX_SHADER
layout(location = 1) in highp vec4 amplitude[NUM];

uniform sampler1D textureData;

float Ampl(int i){
  i = i%DIR_NUM;
  return amplitude[i/4][i%4];
}
#endif

// Define Ampl for fragment shader
#if FRAGMENT_SHADER
in highp vec4 ampl[NUM];

uniform sampler1D textureData;

float Ampl(int i){
  i = i%DIR_NUM;
  return ampl[i/4][i%4];

}
#endif

uniform float profilePeriod;

const float tau = 6.28318530718;

float iAmpl( float angle/*range [0,2pi]*/){
  float a = DIR_NUM*angle/tau;
  int ia = int(floor(a));
  float w = a-ia;
  return (1-w)*Ampl(ia)+w*Ampl((ia+1)%DIR_NUM);
}

vec3 wavePosition(vec3 p){

  vec3 result = vec3(0.0,0.0,0.0);

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kdir_x = dot(p.xy,kdir)+tau*sin(54*a);
    float w = kdir_x/profilePeriod;

    vec4 tt = dx*iAmpl(angle)*texture(textureData,w);

    result.xy += kdir*tt.x;
    result.z += tt.y;
  }

  return result;
}

vec3 waveNormal(vec3 p){

  vec3 tx = vec3(1.0,0.0,0.0);
  vec3 ty = vec3(0.0,1.0,0.0);  

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kdir_x = dot(p.xy,kdir)+tau*sin(54*a);
    float w = kdir_x/profilePeriod;

    vec4 tt = dx*iAmpl(angle)*texture(textureData,w);

    tx.xz += kdir.x*tt.zw;
    ty.yz += kdir.y*tt.zw;
  }

  return normalize(cross(tx,ty));
}

/* 
 * This function is from: https://github.com/hughsk/glsl-hsv2rgb
 */
vec3 hsv2rgb(vec3 c)  {
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 waveColor(vec3 p){

  vec3 color = vec3(0,0,0);
  float norm = 0.001;

  for(int i=0;i<DIR_NUM;i++){
    float angle = (tau*i)/DIR_NUM;
    color += hsv2rgb(vec3(angle,1,Ampl(i)));
    norm += Ampl(i);
  }
  return color/norm;
}
