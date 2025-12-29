#pragma once

#include <image.hpp>
#include <point.hpp>

namespace improcessing {

enum class PerspectiveType {
    kTwoPoint,
    kThreePoint
};

struct SceneObject {
    Point3 center;
    Point3 size;
    Point3 rotation;
    Pixel color;
    PerspectiveType type;
};

} // namespace improcessing