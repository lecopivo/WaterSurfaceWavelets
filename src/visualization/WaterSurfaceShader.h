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
  typedef Attribute<0, Vector3> Position;

  template <int I> using Amplitude = Attribute<1 + I, Vector4>;

  explicit WaterSurfaceShader();

  WaterSurfaceShader &setLightPosition(const Vector3 &position) {
    setUniform(_lightPositionUniform, position);
    return *this;
  }

  WaterSurfaceShader &setAmbientColor(const Color4 &color) {
    setUniform(_ambientColorUniform, color);
    return *this;
  }

  WaterSurfaceShader &setColor(const Color4 &color) {
    setUniform(_colorUniform, color);
    return *this;
  }

  WaterSurfaceShader &setTransformationMatrix(const Matrix4 &matrix) {
    setUniform(_transformationMatrixUniform, matrix);
    return *this;
  }

  WaterSurfaceShader &setNormalMatrix(const Matrix3x3 &matrix) {
    setUniform(_normalMatrixUniform, matrix);
    return *this;
  }

  WaterSurfaceShader &setProjectionMatrix(const Matrix4 &matrix) {
    setUniform(_projectionMatrixUniform, matrix);
    return *this;
  }

  WaterSurfaceShader &setTime(const Float &time) {
    setUniform(_timeUniform, time);
    return *this;
  }

  WaterSurfaceShader &setWaveNumber(const Float &waveNumber) {
    setUniform(_waveNumberUniform, waveNumber);
    return *this;
  }



private:
  Int _lightPositionUniform, _ambientColorUniform, _colorUniform,
      _transformationMatrixUniform, _normalMatrixUniform,
    _projectionMatrixUniform,_timeUniform,_waveNumberUniform;
};

} // namespace Shaders
} // namespace Magnum
