const int NUM_INTEGRATION_NODES = 16*DIR_NUM;

uniform sampler1D profileData;
uniform sampler3D amplitudeData;
uniform float profilePeriod;
uniform float domainSize;
uniform int directionNumber;

const float tau = 6.28318530718;

vec2 positionUV(vec3 pos){
  pos += vec3(domainSize,domainSize,0);
  pos *= 1/(2*domainSize);
  return vec2(pos.x,pos.y);
}

int seed = 40234324;

vec3 wavePosition(vec3 p){

  vec2 uv = positionUV(p);
  vec3 result = vec3(0.0,0.0,0.0);

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kdir_x = dot(p.xy,kdir)+tau*sin(seed*a);
    float w = kdir_x/profilePeriod;

    vec3 pos3 = vec3(uv.x,uv.y,a);
    vec4 tt = dx*texture(amplitudeData,pos3).r*texture(profileData,w);

    result.xy += kdir*tt.x;
    result.z += tt.y;
  }

  return result;
}

vec3 waveNormal(vec3 p){

  vec2 uv = positionUV(p);
  vec3 tx = vec3(1.0,0.0,0.0);
  vec3 ty = vec3(0.0,1.0,0.0);  

  const int N = NUM_INTEGRATION_NODES;
  float da = 1.0/N;
  float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    float angle = a*tau;
    vec2 kdir = vec2(cos(angle),sin(angle));
    float kdir_x = dot(p.xy,kdir)+tau*sin(seed*a);
    float w = kdir_x/profilePeriod;

    vec3 pos3 = vec3(uv.x,uv.y,a);
    vec4 tt = dx*texture(amplitudeData,pos3).r*texture(profileData,w);

    tx.xz += kdir.x*tt.zw;
    ty.yz += kdir.y*tt.zw;
  }

  return normalize(cross(tx,ty));
}
