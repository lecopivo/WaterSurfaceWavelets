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

  //////////////////////////////////////
  lowp vec4 finalAmbientColor = ambientColor;
  //lowp vec4 finalDiffuseColor = vec4(Ampl(0),Ampl(1),Ampl(2),1.0);
  lowp vec4 finalDiffuseColor = color;
  lowp vec4 finalSpecularColor = vec4(1.0,1.0,1.0,1.0);
  lowp vec4 lightColor = vec4(1.0,1.0,1.0,1.0);

  //finalDiffuseColor.rgb += waveColor(pos);

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
