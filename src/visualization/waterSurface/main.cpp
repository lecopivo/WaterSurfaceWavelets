#include <Corrade/Utility/Debug.h>
#include <Corrade/Utility/Resource.h>

using namespace Corrade;

int main() {
    /* Initialize resource manager with name of the group we are using */
    Utility::Resource rs{"data"};

    /* Print out the license */
    Utility::Debug{} << rs.get("WaterSurfaceShader.glsl");
}
