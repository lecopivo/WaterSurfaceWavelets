#pragma once

#include "Magnum/AbstractShaderProgram.h"
#include "Magnum/Math/Color.h"
#include "Magnum/Math/Matrix4.h"
#include "Magnum/Shaders/Generic.h"
#include "Magnum/Shaders/visibility.h"

#include "DirectionNumber.h"

namespace Magnum {
namespace Shaders {

class WaterSurfaceShader : public AbstractShaderProgram {
public:
  typedef Attribute<0, Vector3>        Position;
  
  template<int I>
  using Amplitude = Attribute<1+I, Vector4>;

  explicit WaterSurfaceShader();

  WaterSurfaceShader &setLightPosition(const Vector3 &position);

  WaterSurfaceShader &setAmbientColor(const Color4 &color);

  WaterSurfaceShader &setColor(const Color4 &color);

  WaterSurfaceShader &setTransformationMatrix(const Matrix4 &matrix);

  WaterSurfaceShader &setNormalMatrix(const Matrix3x3 &matrix);

  WaterSurfaceShader &setProjectionMatrix(const Matrix4 &matrix);

private:
  Int _lightPositionUniform, _ambientColorUniform, _colorUniform,
      _transformationMatrixUniform, _normalMatrixUniform,
      _projectionMatrixUniform;
};

} // namespace Shaders
} // namespace Magnum


