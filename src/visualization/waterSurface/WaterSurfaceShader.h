#pragma once

#include "Magnum/AbstractShaderProgram.h"
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

  WaterSurfaceShader &setProfilePeriod(const Float &profilePeriod) {
    setUniform(_profilePeriodUniform, profilePeriod);
    return *this;
  }

  WaterSurfaceShader &setDomainSize(const Float &size) {
    setUniform(_domainSizeUniform, size);
    return *this;
  }

  WaterSurfaceShader &setDirectionNumber(const Int &dir_num) {
    setUniform(_directionNumberUniform, dir_num);
    return *this;
  }

  WaterSurfaceShader &bindProfileTexture(Texture1D &texture) {
    texture.bind(ProfileTextureLayer);
    return *this;
  }

  WaterSurfaceShader &bindAmplitudeTexture(Texture3D &texture) {
    texture.bind(AmplitudeTextureLayer);
    return *this;
  }

private:
  enum : Int { ProfileTextureLayer = 0, AmplitudeTextureLayer = 1 };

  Int _lightPositionUniform, _ambientColorUniform, _colorUniform,
      _transformationMatrixUniform, _normalMatrixUniform,
      _projectionMatrixUniform, _profilePeriodUniform, _domainSizeUniform,
      _directionNumberUniform;
};

} // namespace Shaders
} // namespace Magnum
