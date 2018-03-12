/* DIR_NUM is defined in DirectionNumber.h which is included during shader compilation because #include is not supported for some fucking stupid reson@&#! */
const int NUM = DIR_NUM/4;

uniform lowp vec4 ambientColor;
uniform lowp vec4 color;
uniform mat3 normalMatrix;

in mediump vec3 transformedNormal;
in highp vec3 lightDirection;
in highp vec3 cameraDirection;
in highp vec3 pos;
in highp vec4 ampl[NUM];

layout(location = 0) out lowp vec4 fragmentColor;

float Ampl(int i){
  return ampl[i/4][i%4];
}

const float id2rad = 6.28318530718/DIR_NUM;

float sine_wave(float a, float k, float kx, float t) {

  const float g = 9.81;

  float omega = sqrt(g / k);
  float theta = kx - omega * t;
  float nu    = a * cos(theta);
  
  return nu;
}

float sine_wave_der(float a, float k, float kx, float t) {

  const float g = 9.81;

  float omega = sqrt(g / k);
  float theta = kx - omega * t;
  float der    = - a * sin(theta);
  
  return der;
}


void main() {
  highp vec3 p = pos;
  highp vec3 normal;

  highp vec3 tx = vec3(1.0,0.0,0.0);
  highp vec3 ty = vec3(0.0,1.0,0.0);
  for(int i=0;i<DIR_NUM;i++){
    float k = 1;
    float angle = i*id2rad;
    vec2 kvec = k*vec2(cos(angle),sin(angle));
    float kx = dot(p.xy,kvec);
    
    p.z += sine_wave(Ampl(i),k,kx,0);
    
    float der = sine_wave_der(Ampl(i),k,kx,0)+sine_wave_der(0.05*Ampl(i),100*k,100*kx,0);
    tx.z += kvec.x*der;
    ty.z += kvec.y*der;
  }
  
  normal = normalize(normalMatrix*cross(tx,ty));

  //////////////////////////////////////
  lowp vec4 finalAmbientColor = ambientColor;
  lowp vec4 finalDiffuseColor = color;
  lowp vec4 finalSpecularColor = vec4(1.0,1.0,1.0,1.0);
  lowp vec4 lightColor = vec4(1.0,1.0,1.0,1.0);

  /* Ambient color */
  fragmentColor = finalAmbientColor;

  mediump vec3 normalizedTransformedNormal = normalize(normal);
  highp vec3 normalizedLightDirection = normalize(lightDirection);

  /* Add diffuse color */
  lowp float intensity = max(0.0, dot(normalizedTransformedNormal, normalizedLightDirection));
  fragmentColor += finalDiffuseColor*lightColor*intensity;

  /* Add specular color, if needed */
  if(intensity > 0.001) {
    highp vec3 reflection = reflect(-normalizedLightDirection, normalizedTransformedNormal);
    highp float shininess = 80.0;
    mediump float specularity = pow(max(0.0, dot(normalize(cameraDirection), reflection)), shininess);
    fragmentColor += finalSpecularColor*specularity;
  }

  //fragmentColor = vec4(1.0,0.0,0.0,1.0);
}
