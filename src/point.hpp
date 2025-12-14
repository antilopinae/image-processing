#pragma once

#include <cstddef>
#include <cmath>
#include <iomanip>
#include <cassert>

namespace improcessing {
    template<typename T>
    struct PointBase {
        T x, y;
    };

    using Point2D = PointBase<size_t>;

    struct Point : PointBase<double> {
        constexpr static double EPS = 1e-9;

        Point() : PointBase(0.0, 0.0) {
        }

        Point(double x, double y) : PointBase(x, y) {
        }

        [[nodiscard]] auto Dot(const Point &o) const noexcept -> double { return x * o.x + y * o.y; }
        [[nodiscard]] auto Cross(const Point &o) const noexcept -> double { return x * o.y - y * o.x; }
        [[nodiscard]] auto Len() const noexcept -> double { return hypot(x, y); }

        [[nodiscard]] auto Normalized() const noexcept -> Point {
            double len = Len();
            if (len <= EPS) return {0, 0};
            return {x / len, y / len};
        }

        [[nodiscard]] auto Equal(const Point &o, double tol = 1e-9) const noexcept -> bool {
            return fabs(x - o.x) <= tol && fabs(y - o.y) <= tol;
        }

        [[nodiscard]] auto ToPoint2D() const -> Point2D {
            assert(x >= 0 && y >= 0 && "Invalid Point2D");
            return Point2D{
                static_cast<size_t>(std::floor(x)),
                static_cast<size_t>(std::floor(y))
            };
        }
    };

    inline auto operator+(const Point &a, const Point &b) noexcept -> Point {
        return {a.x + b.x, a.y + b.y};
    }

    inline auto operator-(const Point &a, const Point &b) noexcept -> Point {
        return {a.x - b.x, a.y - b.y};
    }

    inline auto operator*(const Point &p, double k) noexcept -> Point {
        return {p.x * k, p.y * k};
    }

    inline auto operator*(double k, const Point &p) noexcept -> Point {
        return {p.x * k, p.y * k};
    }

    inline auto operator/(const Point &p, double k) noexcept -> Point {
        return {p.x / k, p.y / k};
    }

    inline auto operator<<(std::ostream &os, const Point &p) -> std::ostream & {
        os.setf(std::ios::fixed);
        os << std::setprecision(6) << "(" << p.x << ", " << p.y << ")";
        os.unsetf(std::ios::fixed);
        return os;
    }

    /*!
     * @brief Classification of a point relative to an oriented segment a->b
     * @param a, b Coordinates of oriented segment
     * @param p Coordinate of point
     * @return +1 if the point is on the left, -1 if on the RIGHT, 0 if on the line
     */
    inline auto ClassifyPointToEdge(const Point &a, const Point &b, const Point &p) -> int {
        double cr = (b - a).Cross(p - a);
        if (cr > Point::EPS) return +1;
        if (cr < -Point::EPS) return -1;
        return 0;
    }

    struct Point3 : Point {
        double z;

        Point3() : Point(0, 0), z(0) {
        }

        Point3(double x, double y, double z) : Point(x, y), z(z) {
        }

        [[nodiscard]] auto Dot(const Point3 &o) const noexcept -> double { return x * o.x + y * o.y + z * o.z; }

        [[nodiscard]] auto Cross(const Point3 &o) const noexcept -> Point3 {
            return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
        }

        [[nodiscard]] auto Len() const noexcept -> double { return std::sqrt(x * x + y * y + z * z); }

        [[nodiscard]] auto Normalized() const noexcept -> Point3 {
            double len = Len();
            if (len <= EPS) return {0, 0, 0};
            return {x / len, y / len, z / len};
        }
    };

    inline auto operator+(const Point3 &a, const Point3 &b) noexcept -> Point3 {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    }

    inline auto operator-(const Point3 &a, const Point3 &b) noexcept -> Point3 {
        return {a.x - b.x, a.y - b.y, a.z - b.z};
    }

    inline auto operator*(const Point3 &p, double k) noexcept -> Point3 {
        return {p.x * k, p.y * k, p.z * k};
    }

    inline auto operator*(double k, const Point3 &p) noexcept -> Point3 {
        return p * k;
    }

    inline auto operator/(const Point3 &p, double k) noexcept -> Point3 {
        return {p.x / k, p.y / k, p.z / k};
    }
} // namespace improcessing
