/* DIR_NUM is defined in DirectionNumber.h which is included during shader compilation because #include is not supported for some fucking stupid reson@&#! */
const int NUM = DIR_NUM/4;

uniform lowp vec4 ambientColor;
uniform lowp vec4 color;

in mediump vec3 transformedNormal;
in highp vec3 lightDirection;
in highp vec3 cameraDirection;
in highp vec4 pos;
in highp vec4 ampl[NUM];

layout(location = 0) out lowp vec4 fragmentColor;

float Ampl(int i){
  return ampl[i/4][i%4];
}

void main() {
  lowp vec4 finalAmbientColor = ambientColor;
  lowp vec4 finalDiffuseColor = color;
  lowp vec4 finalSpecularColor = vec4(1.0,1.0,1.0,1.0);
  lowp vec4 lightColor = vec4(1.0,1.0,1.0,1.0);

  /* Ambient color */
  fragmentColor = finalAmbientColor;

  mediump vec3 normalizedTransformedNormal = normalize(transformedNormal);
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

  //fragmentColor = vec4(ampl[0],pos[1],pos[2],1.0);
  //fragmentColor = vec4(1.0,0.0,0.0,1.0);
}
