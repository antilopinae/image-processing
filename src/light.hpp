#pragma once

#include <image.hpp>
#include <point.hpp>

namespace improcessing {

struct Light {
    Point3 position;
    double ambient, diffuse, specular;
    Pixel color;
};

} // namespace improcessing