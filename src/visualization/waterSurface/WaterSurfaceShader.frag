// The final shader is in several files and the final version is stitched together in the constructor of the class WaterSurfaceShader. This is because GLSL does not provide #include directive
//
// The final fragmetn shader contain:
// #include "DirectionNumber.h"
// # define FRAGMENT_SHADER
// #include "WaterSurfaceShader.glsl"

uniform lowp vec4 ambientColor;
uniform lowp vec4 color;
uniform mat3 normalMatrix;
uniform float time;
uniform float waveNumber;

in mediump vec3 transformedNormal;
in highp vec3 lightDirection;
in highp vec3 cameraDirection;
in highp vec3 pos;

layout(location = 0) out lowp vec4 fragmentColor;

void main() {
  highp vec3 normal = normalMatrix*waveNormal(pos);
  //highp vec3 normal = normalMatrix*vec3(0,0,1);

  //////////////////////////////////////
  lowp vec4 finalAmbientColor = ambientColor;
  lowp vec4 finalDiffuseColor = color;
  lowp vec4 finalSpecularColor = vec4(1.0,1.0,1.0,1.0);
  lowp vec4 lightColor = vec4(1.0,1.0,1.0,1.0);

  if(pos.x < -50 || pos.x > 50 || pos.y<-50 || pos.y>50)
    finalDiffuseColor.rgb = vec3(0.6,0.6,0.6);

  /* Ambient color */
  fragmentColor = finalAmbientColor;

  mediump vec3 normalizedTransformedNormal = normalize(normal);
  highp vec3 normalizedLightDirection = normalize(lightDirection);

  /* Add diffuse color */
  lowp float intensity = max(0.0, dot(normalizedTransformedNormal, normalizedLightDirection));
  fragmentColor += finalDiffuseColor*lightColor*intensity;

  /* Add specular color, if needed */
  if(intensity > 0.001) {
    // vec3 ref = reflect(normalize(cameraDirection),normalizedTransformedNormal);
    // highp float sky = max(0,pow((1-abs(dot(normalize(cameraDirection),normalizedTransformedNormal))),1)*sin(20*ref.x)*sin(20*ref.y)*sin(20*ref.z));
    highp vec3 reflection = reflect(-normalizedLightDirection, normalizedTransformedNormal);
    highp float shininess = 80.0;
    mediump float specularity = pow(max(0.0, dot(normalize(cameraDirection), reflection)), shininess);
    fragmentColor += finalSpecularColor*specularity;// + sky*vec4(1,1,1,1);
  }
}
