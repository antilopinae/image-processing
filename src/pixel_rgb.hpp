#pragma once

#include <point.hpp>

namespace improcessing {
template<typename T>
struct PixelRGB {
    union {
        struct {
            T r;
            T g;
            T b;
        };

        T data[3];
    };

    explicit PixelRGB() : r(0), g(0), b(0) {}

    PixelRGB(T r_, T g_, T b_) : r(r_), g(g_), b(b_) {}

    bool operator==(const PixelRGB& other) const
    {
        return r == other.r && g == other.g && b == other.b;
    }

    bool operator!=(const PixelRGB& other) const
    {
        return !(*this == other);
    }

    template<typename Scalar>
    PixelRGB operator*(Scalar scalar) const
    {
        return PixelRGB(static_cast<T>(r * scalar), static_cast<T>(g * scalar), static_cast<T>(b * scalar));
    }

    template<typename Scalar>
    PixelRGB& operator*=(Scalar scalar)
    {
        r = static_cast<T>(r * scalar);
        g = static_cast<T>(g * scalar);
        b = static_cast<T>(b * scalar);
        return *this;
    }

    PixelRGB operator+(const PixelRGB& other) const
    {
        return PixelRGB(r + other.r, g + other.g, b + other.b);
    }

    PixelRGB operator*(const PixelRGB& other) const
    {
        return PixelRGB(r * other.r, g * other.g, b * other.b);
    }

    Point3 ToPoint3() const
    {
        return Point3(static_cast<double>(r) / 255.0, static_cast<double>(g) / 255.0, static_cast<double>(b) / 255.0);
    }

    static PixelRGB<unsigned char> FromPoint3(const Point3& p)
    {
        return PixelRGB<unsigned char>(static_cast<unsigned char>(std::min(1.0, std::max(0.0, p.x)) * 255.0),
                                       static_cast<unsigned char>(std::min(1.0, std::max(0.0, p.y)) * 255.0),
                                       static_cast<unsigned char>(std::min(1.0, std::max(0.0, p.z)) * 255.0));
    }
};

template<typename T, typename Scalar>
PixelRGB<T> operator*(Scalar scalar, const PixelRGB<T>& pixel)
{
    return pixel * scalar;
}
} // namespace improcessing
