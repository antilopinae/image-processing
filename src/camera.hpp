#pragma once

#include <point.hpp>

namespace improcessing {

class Camera {
public:
    Point3 position{0, 0, -10};
    double focal_length = 500.0;
    double near_plane   = 0.1;

    Point2D Project(const Point3& p, size_t width, size_t height) const
    {
        double z_eff = p.z - position.z;

        if (z_eff < 0.1) {
            z_eff = 0.1;
        }

        double x_proj = (p.x * focal_length) / z_eff;
        double y_proj = (p.y * focal_length) / z_eff;

        return {static_cast<size_t>((width / 2.0) + x_proj), static_cast<size_t>((height / 2.0) - y_proj)};
    }
};

} // namespace improcessing