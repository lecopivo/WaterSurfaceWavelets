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


highp float iAmpl(highp float a/*range [0,1]*/){
  a *= DIR_NUM;
  highp int ia = int(floor(a));
  highp float w = a-ia;
  return (1-w)*Ampl(ia)+w*Ampl((ia+1)%DIR_NUM);
}

const highp float tau = 6.28318530718;

highp vec2 sine_wave(highp float a,highp float k,highp float kx,highp float t) {

  const highp float g = 9.81;
  
  highp float omega = sqrt(g / k);
  highp float theta = kx - omega * t;
  highp vec2 nu = a*vec2(cos(theta),-k*sin(theta));
  
  return nu;
}

// returns (normal.x,normal.y,height)
// actual normal is normal = normalize(vec3(normal.x,normal.y,1))
highp vec3 wave(highp vec3 p, highp float k, highp float t){

  highp vec3 result = vec3(0.0,0.0,0.0);

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kx = k*dot(p.xy,kdir);

    vec2 sw = dx*(sine_wave(iAmpl(a),k,kx+tau*sin(54*a),t)+sine_wave(0.05*iAmpl(a),5*k,5*kx,t)+sine_wave(0.01*iAmpl(a),10*k,10*kx,t));

    result.xy -= kdir*sw.y;
    result.z += sw.x;
  }

  result.z = 10*texture(textureData,0.1*p.x).y;

  return result;
}


vec3 wavePosition(vec3 p, float k){

  vec3 result = vec3(0.0,0.0,0.0);

  float lambda = tau/k;

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kx = k*dot(p.xy,kdir)+tau*sin(54*a);
    float w = kx/lambda;

    vec4 tt = dx*iAmpl(a)*texture(textureData,w);

    result.xy += kdir*tt.x;
    result.z += tt.y;
  }

  return result;
}

vec3 waveNormal(vec3 p, float k){

  vec3 tx = vec3(1.0,0.0,0.0);
  vec3 ty = vec3(0.0,1.0,0.0);  

  float lambda = tau/k;

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kx = k*dot(p.xy,kdir)+tau*sin(54*a);
    float w = kx/lambda;

    vec4 tt = dx*iAmpl(a)*texture(textureData,w);

    tx.xz += kdir.x*tt.zw;
    ty.yz += kdir.y*tt.zw;
  }

  return normalize(cross(tx,ty));
}
