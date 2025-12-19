#pragma once

namespace improcessing {
    struct Matrix4x4 {
        double m[4][4] = {0};

        static Matrix4x4 Identity() {
            Matrix4x4 res;
            for (int i = 0; i < 4; ++i) res.m[i][i] = 1.0;
            return res;
        }

        static Matrix4x4 Rotation(Point3 axis, double angle) {
            axis = axis.Normalized();
            double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
            Matrix4x4 res = Identity();
            res.m[0][0] = t * axis.x * axis.x + c;
            res.m[0][1] = t * axis.x * axis.y - s * axis.z;
            res.m[0][2] = t * axis.x * axis.z + s * axis.y;
            res.m[1][0] = t * axis.x * axis.y + s * axis.z;
            res.m[1][1] = t * axis.y * axis.y + c;
            res.m[1][2] = t * axis.y * axis.z - s * axis.x;
            res.m[2][0] = t * axis.x * axis.z - s * axis.y;
            res.m[2][1] = t * axis.y * axis.z + s * axis.x;
            res.m[2][2] = t * axis.z * axis.z + c;
            return res;
        }

        Point3 Transform(const Point3 &p) const {
            double x = p.x * m[0][0] + p.y * m[0][1] + p.z * m[0][2] + m[0][3];
            double y = p.x * m[1][0] + p.y * m[1][1] + p.z * m[1][2] + m[1][3];
            double z = p.x * m[2][0] + p.y * m[2][1] + p.z * m[2][2] + m[2][3];
            double w = p.x * m[3][0] + p.y * m[3][1] + p.z * m[3][2] + m[3][3];
            if (std::abs(w) > 1e-9) return {x / w, y / w, z / w};
            return {x, y, z};
        }
    };
}
