#pragma once

#include "Magnum/GL/AbstractShaderProgram.h"
#include "Magnum/Math/Color.h"
#include "Magnum/Math/Matrix4.h"
#include "Magnum/Shaders/Generic.h"
#include "Magnum/Shaders/visibility.h"
#include "Magnum/Texture.h"

#include "DirectionNumber.h"

namespace Magnum {
namespace Shaders {

class WaterSurfaceShader : public AbstractShaderProgram {
public:
  typedef GL::Attribute<0, Vector3> Position;

  template <int I> using Amplitude = GL::Attribute<1 + I, Vector4>;

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

  WaterSurfaceShader &setProfilePeriod(const Float &profilePeriod) {
    setUniform(_profilePeriodUniform, profilePeriod);
    return *this;
  }

  WaterSurfaceShader &setGerstnerParameter(const Float &gerstnerParameter) {
    setUniform(_gerstnerParameterUniform, gerstnerParameter);
    return *this;
  }

  WaterSurfaceShader &setWaveDirectionToShow(const int &direction) {
    setUniform(_waveDirectionToShowUniform, direction);
    return *this;
  }

  WaterSurfaceShader &bindTexture(GL::Texture1D &texture) {
    texture.bind(TextureLayer);
    return *this;
  }

private:
  enum : Int { TextureLayer = 0 };

  Int _lightPositionUniform, _ambientColorUniform, _colorUniform,
      _transformationMatrixUniform, _normalMatrixUniform,
      _projectionMatrixUniform, _timeUniform, _profilePeriodUniform,
      _gerstnerParameterUniform, _waveDirectionToShowUniform;
};

} // namespace Shaders
} // namespace Magnum
