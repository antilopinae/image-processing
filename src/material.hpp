#pragma once

namespace improcessing {

struct Material {
    double k_ambient, k_diffuse, k_specular;
    double shininess;
};

} // namespace improcessing