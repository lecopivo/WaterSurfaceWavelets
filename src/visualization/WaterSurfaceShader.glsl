const int NUM = DIR_NUM/4;

// Define Ampl for vertex shader
#if VERTEX_SHADER
layout(location = 1) in highp vec4 amplitude[NUM];

float Ampl(int i){
  i = i%DIR_NUM;
  return amplitude[i/4][i%4];
}
#endif

// Define Ampl for fragment shader
#if FRAGMENT_SHADER
in highp vec4 ampl[NUM];

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

  int N = 32;
  highp float da = 1.0/N;
  highp float dx = DIR_NUM*tau/N;
  for(float a = 0;a<1;a+=da){

    highp float angle = a*tau;
    highp vec2 kdir = vec2(cos(angle),sin(angle));
    highp float kx = k*dot(p.xy,kdir);

    highp vec2 sw = dx*(sine_wave(iAmpl(a),k,kx+tau*sin(54*a),t)+sine_wave(0.05*iAmpl(a),5*k,5*kx,t)+sine_wave(0.01*iAmpl(a),10*k,10*kx,t));

    result.xy -= kdir*sw.y;
    result.z += sw.x;
  }

  return result;
}
