uniform highp mat4 transformationMatrix;
uniform highp mat4 projectionMatrix;
uniform mediump mat3 normalMatrix;
uniform highp vec3 light;
uniform float time;

layout(location = 0) in highp vec4 position;

out mediump vec3 transformedNormal;
out highp vec3 lightDirection;
out highp vec3 cameraDirection;
out highp vec3 pos;
out highp vec4 ampl[NUM];

void main() {
  highp vec3 p = position.xyz/position.w;

  p += wavePosition(p);
  vec3 normal = vec3(0.0,0.0,1.0);
  
  /* Transformed vertex position */
  highp vec4 transformedPosition4 = transformationMatrix*vec4(p.x,p.y,p.z,1);
  highp vec3 transformedPosition = transformedPosition4.xyz/transformedPosition4.w;

  /* Transformed normal vector */
  transformedNormal = normalMatrix*normal;

  /* Direction to the light */
  lightDirection = normalize(light - transformedPosition);

  /* Direction to the camera */
  cameraDirection = -transformedPosition;

  /* Transform the position */
  gl_Position = projectionMatrix*transformedPosition4;

  pos = position.xyz/position.w;

  for(int i=0;i<NUM;i++){
    ampl[i] = amplitude[i];
  }
}
