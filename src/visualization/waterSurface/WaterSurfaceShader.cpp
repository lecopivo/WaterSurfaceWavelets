#include "WaterSurfaceShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Resource.h>

#include <Magnum/Extensions.h>
#include <Magnum/Shader.h>
#include <Magnum/Texture.h>

namespace Magnum {
namespace Shaders {

WaterSurfaceShader::WaterSurfaceShader() {

  Utility::Resource rs("WaterSurface");

  GL::Shader vert{Version::GL330, Shader::Type::Vertex},
      frag{Version::GL330, Shader::Type::Fragment};

  vert.addSource(rs.get("DirectionNumber.h"))
      .addSource("#define VERTEX_SHADER 1\n")
      .addSource(rs.get("WaterSurfaceShader.glsl"))
      .addSource(rs.get("WaterSurfaceShader.vert"));
  frag.addSource(rs.get("DirectionNumber.h"))
      .addSource("#define FRAGMENT_SHADER 1\n")
      .addSource(rs.get("WaterSurfaceShader.glsl"))
      .addSource(rs.get("WaterSurfaceShader.frag"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, frag}));
  attachShaders({vert, frag});
  CORRADE_INTERNAL_ASSERT(link());

  _lightPositionUniform        = uniformLocation("light");
  _ambientColorUniform         = uniformLocation("ambientColor");
  _colorUniform                = uniformLocation("color");
  _transformationMatrixUniform = uniformLocation("transformationMatrix");
  _projectionMatrixUniform     = uniformLocation("projectionMatrix");
  _normalMatrixUniform         = uniformLocation("normalMatrix");
  _timeUniform                 = uniformLocation("time");
  _profilePeriodUniform        = uniformLocation("profilePeriod");
  _gerstnerParameterUniform    = uniformLocation("gerstnerParameter");
  _waveDirectionToShowUniform  = uniformLocation("waveDirectionToShow");

  setUniform(uniformLocation("textureData"), TextureLayer);
}

} // namespace Shaders
} // namespace Magnum
