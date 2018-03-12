#include "WaterSurfaceShader.h"

#include <Corrade/Utility/Resource.h>

//#include "Magnum/Context.h"
#include "Magnum/Extensions.h"
#include "Magnum/Shader.h"
#include "Magnum/Texture.h"

namespace Magnum {
namespace Shaders {

WaterSurfaceShader &
WaterSurfaceShader::setLightPosition(const Vector3 &position) {
  setUniform(_lightPositionUniform, position);
  return *this;
}

WaterSurfaceShader &WaterSurfaceShader::setAmbientColor(const Color4 &color) {
  setUniform(_ambientColorUniform, color);
  return *this;
}

WaterSurfaceShader &WaterSurfaceShader::setColor(const Color4 &color) {
  setUniform(_colorUniform, color);
  return *this;
}

WaterSurfaceShader &
WaterSurfaceShader::setTransformationMatrix(const Matrix4 &matrix) {
  setUniform(_transformationMatrixUniform, matrix);
  return *this;
}

WaterSurfaceShader &
WaterSurfaceShader::setNormalMatrix(const Matrix3x3 &matrix) {
  setUniform(_normalMatrixUniform, matrix);
  return *this;
}

WaterSurfaceShader &
WaterSurfaceShader::setProjectionMatrix(const Matrix4 &matrix) {
  setUniform(_projectionMatrixUniform, matrix);
  return *this;
}

WaterSurfaceShader::WaterSurfaceShader() {

  Utility::Resource rs("WaterSurfaceShaders");

  Shader vert{Version::GL330, Shader::Type::Vertex},
      frag{Version::GL330, Shader::Type::Fragment};

  vert.addSource(rs.get("DirectionNumber.h"))
      .addSource(rs.get("WaterSurfaceShader.vert"));
  frag.addSource(rs.get("DirectionNumber.h"))
      .addSource(rs.get("WaterSurfaceShader.frag"));

  CORRADE_INTERNAL_ASSERT(Shader::compile({vert, frag}));
  attachShaders({vert, frag});
  CORRADE_INTERNAL_ASSERT(link());

  _lightPositionUniform        = uniformLocation("light");
  _ambientColorUniform         = uniformLocation("ambientColor");
  _colorUniform                = uniformLocation("color");
  _transformationMatrixUniform = uniformLocation("transformationMatrix");
  _projectionMatrixUniform     = uniformLocation("projectionMatrix");
  _normalMatrixUniform         = uniformLocation("normalMatrix");
}

} // namespace Shaders
} // namespace Magnum
