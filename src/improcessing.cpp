#include <cmath>
#include <cstdlib>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <image.hpp>
#include <improcessing.hpp>
#include <matrix4x4.hpp>

#include <file_closer.hpp>
#include <png_guard.hpp>

namespace improcessing {
auto ReadImage(const std::string &filename) -> std::expected<Image, boost::system::error_code>
{
    std::unique_ptr<FILE, details::FileCloser> file(std::fopen(filename.c_str(), "rb"));
    if (!file) {
        return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)};
    }

    PngGuard png(false);

    png_init_io(png.png_ptr, file.get());
    png_read_info(png.png_ptr, png.info_ptr);

    png_uint_32 WidthGray  = 0;
    png_uint_32 HeightGray = 0;
    int bit_depth          = 0;
    int color_type         = 0;

    png_get_IHDR(png.png_ptr, png.info_ptr, &WidthGray, &HeightGray, &bit_depth, &color_type, NULL, NULL, NULL);

    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png.png_ptr);
    }

    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
        png_set_rgb_to_gray_fixed(png.png_ptr, 1, -1, -1);
    }

    if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_strip_alpha(png.png_ptr);
    }

    if (bit_depth == 16) {
        png_set_strip_16(png.png_ptr);
    }

    png_read_update_info(png.png_ptr, png.info_ptr);

    Image image(WidthGray, HeightGray);

    auto row_ptrs = image.GetGrayPtrs();

    png_read_image(png.png_ptr, row_ptrs.data());
    png_read_end(png.png_ptr, png.info_ptr);

    return std::move(image);
}

auto SaveImage(const std::string &filename, const Image &image, Image::Type type)
    -> std::expected<void, boost::system::error_code>
{
    std::unique_ptr<FILE, details::FileCloser> file(std::fopen(filename.c_str(), "wb"));
    if (!file) {
        return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::io_error)};
    }

    PngGuard png(true);

    png_init_io(png.png_ptr, file.get());

    switch (type) {
        case Image::Type::kGray: {
            png_set_IHDR(png.png_ptr,
                         png.info_ptr,
                         image.WidthGray(),
                         image.HeightGray(),
                         8,
                         PNG_COLOR_TYPE_GRAY,
                         PNG_INTERLACE_NONE,
                         PNG_COMPRESSION_TYPE_DEFAULT,
                         PNG_FILTER_TYPE_DEFAULT);

            png_write_info(png.png_ptr, png.info_ptr);

            auto row_ptrs = image.GetGrayPtrs();

            png_write_image(png.png_ptr, row_ptrs.data());
        } break;
        case Image::Type::kRGB: {
            png_set_IHDR(png.png_ptr,
                         png.info_ptr,
                         image.WidthRgb(),
                         image.HeightRgb(),
                         8,
                         PNG_COLOR_TYPE_RGB,
                         PNG_INTERLACE_NONE,
                         PNG_COMPRESSION_TYPE_DEFAULT,
                         PNG_FILTER_TYPE_DEFAULT);

            png_write_info(png.png_ptr, png.info_ptr);

            auto row_ptrs = image.GetGrayPtrs();

            png_write_image(png.png_ptr, row_ptrs.data());
        } break;
    }

    png_write_end(png.png_ptr, png.info_ptr);

    return {};
}

auto MakeCircularGrayscale(Image &image, double radius_fraction) -> std::expected<void, boost::system::error_code>
{
    const auto cx             = (image.WidthGray() - 1) / 2.0;
    const auto cy             = (image.HeightGray() - 1) / 2.0;
    const auto r              = std::min(image.WidthGray(), image.HeightGray()) * radius_fraction;
    const auto edge_WidthGray = std::max(1.0, r * 0.02);

    for (int y = 0; y < image.HeightGray(); ++y) {
        for (int x = 0; x < image.WidthGray(); ++x) {
            const auto dx = x - cx;
            const auto dy = y - cy;
            const auto d  = std::sqrt(dx * dx + dy * dy);

            uint8_t value = 0;

            if (d <= r - edge_WidthGray) {
                const auto t = d / (r - edge_WidthGray);
                auto v       = 255.0 * (1.0 - 0.7 * t);

                v     = std::clamp(v, 0.0, 255.0);
                value = static_cast<uint8_t>(std::round(v));
            } else if (d <= r + edge_WidthGray) {
                const auto t = (d - (r - edge_WidthGray)) / (2.0 * edge_WidthGray);
                auto v       = 255.0 * (1.0 - t);

                v     = std::clamp(v, 0.0, 255.0);
                value = static_cast<uint8_t>(std::round(v));
            }

            image.DataGray()[y * image.WidthGray() + x] = value;
        }
    }

    return {};
}

auto Blend(const Image &A, const Image &B, const Image &Alpha) -> std::expected<Image, boost::system::error_code>
{
    if (A.WidthGray() > B.WidthGray() || A.HeightGray() > B.HeightGray() || A.WidthGray() > Alpha.WidthGray() ||
        A.HeightGray() > Alpha.HeightGray()) {
        return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::invalid_argument)};
    }

    Image out;
    out.ResizeGray(A.WidthGray(), A.HeightGray());

    const auto total = A.WidthGray() * A.HeightGray();

    for (size_t i = 0; i < total; ++i) {
        const auto a     = A.DataGray()[i];
        const auto b     = B.DataGray()[i];
        const auto alpha = Alpha.DataGray()[i];

        const auto blended = static_cast<uint64_t>(alpha) * b + static_cast<uint64_t>(255 - alpha) * a;

        out.DataGray()[i] = static_cast<uint8_t>(blended / 255);
    }

    return std::move(out);
}

auto FloydSteinbergDither(Image &img, uint8_t n_levels) -> std::expected<void, boost::system::error_code>
{
    const int rows   = img.HeightGray();
    const int cols   = img.WidthGray();
    const int levels = 1 << n_levels;

    struct {
        Image::gray_type newval;
    } lut[256];

    for (int i = 0; i < 256; ++i) {
        int nearest   = (i * (levels - 1) + 127) / 255;
        int quantized = (nearest * 255 + (levels - 1) / 2) / (levels - 1);
        lut[i]        = {static_cast<Image::gray_type>(quantized)};
    }

    int err;
    auto apply_error = [&](int x, int y, int num, int den) {
        if (y >= 0 && y < rows && x >= 0 && x < cols) {
            int val   = img(x, y) + (err * num) / den;
            img(x, y) = std::clamp(val, 0, 255);
        }
    };

    for (int y = 0; y < rows; ++y) {
        const bool reverse = y & 1;
        const auto start   = reverse ? cols - 1 : 0;
        const auto end     = reverse ? -1 : cols;
        const auto step    = reverse ? -1 : 1;

        for (int x = start; x != end; x += step) {
            const auto old_val = img(x, y);
            const auto entry   = lut[old_val];
            const auto new_val = entry.newval;
            img(x, y)          = new_val;

            err = old_val - new_val;

            if (!reverse) {
                apply_error(x + 1, y, 7, 16);
                apply_error(x - 1, y + 1, 3, 16);
                apply_error(x, y + 1, 5, 16);
                apply_error(x + 1, y + 1, 1, 16);
            } else {
                apply_error(x - 1, y, 7, 16);
                apply_error(x + 1, y + 1, 3, 16);
                apply_error(x, y + 1, 5, 16);
                apply_error(x - 1, y + 1, 1, 16);
            }
        }
    }

    return {};
}

auto DrawLine(Image &img, Point2D start, Point2D end, Pixel color) -> std::expected<void, boost::system::error_code>
{
    int x = start.x, y = start.y;
    int dx = end.x - start.x, dy = end.y - start.y;
    int ix, iy, e, i;
    if (dx > 0) {
        ix = 1;
    } else if (dx < 0) {
        ix = -1;
        dx = -dx;
    } else {
        ix = 0;
    }

    if (dy > 0) {
        iy = 1;
    } else if (dy < 0) {
        iy = -1;
        dy = -dy;
    } else {
        iy = 0;
    }

    if (dx >= dy) {
        e = 2 * dy - dx;
        if (iy >= 0) {
            for (i = 0; i <= dx; i++) {
                img.GetRGBPixel(x, y) = color;
                if (e >= 0) {
                    y += iy;
                    e -= 2 * dx;
                }
                x += ix;
                e += dy * 2;
            }
        } else {
            for (i = 0; i <= dx; i++) {
                img.GetRGBPixel(x, y) = color;
                if (e > 0) {
                    y += iy;
                    e -= 2 * dx;
                }
                x += ix;
                e += dy * 2;
            }
        }
    } else {
        e = 2 * dx - dy;
        if (ix >= 0) {
            for (i = 0; i <= dy; i++) {
                img.GetRGBPixel(x, y) = color;
                if (e >= 0) {
                    x += ix;
                    e -= 2 * dy;
                }
                y += iy;
                e += dx * 2;
            }
        } else {
            for (i = 0; i <= dy; i++) {
                img.GetRGBPixel(x, y) = color;
                if (e > 0) {
                    x += ix;
                    e -= 2 * dy;
                }
                y += iy;
                e += dx * 2;
            }
        }
    }

    return {};
}

auto DrawPolygonEdges(Image &img, const std::vector<Point2D> &poly, Pixel color)
    -> std::expected<void, boost::system::error_code>
{
    auto size = poly.size();

    for (auto i = 0u; i < size; ++i) {
        if (i == size - 1) {
            if (auto res = DrawLine(img, poly[i], poly[0], color); !res) {
                return std::unexpected{res.error()};
            }
            break;
        }

        if (auto res = DrawLine(img, poly[i], poly[i + 1], color); !res) {
            return std::unexpected{res.error()};
        }
    }

    return {};
}

auto FillCircle(Image &img, Point center, double radius, Pixel color) -> void
{
    int r  = static_cast<int>(std::ceil(radius));
    int cx = static_cast<int>(center.x);
    int cy = static_cast<int>(center.y);

    for (int y = -r; y <= r; ++y) {
        for (int x = -r; x <= r; ++x) {
            if (x * x + y * y <= radius * radius) {
                int px = cx + x;
                int py = cy + y;
                if (px >= 0 && px < (int)img.WidthRgb() && py >= 0 && py < (int)img.HeightRgb()) {
                    img.GetRGBPixel(px, py) = color;
                }
            }
        }
    }
}

void DrawThickLine(Image &img, Point2D start2d, Point2D end2d, double thickness, Pixel color, LineCap cap)
{
    Point A(static_cast<double>(start2d.x), static_cast<double>(start2d.y));
    Point B(static_cast<double>(end2d.x), static_cast<double>(end2d.y));

    if (A.Equal(B)) {
        if (cap == LineCap::kRound) {
            FillCircle(img, A, thickness / 2.0, color);
        } else {
            img.GetRGBPixel(start2d.x, start2d.y) = color;
        }
        return;
    }

    Point dir = (B - A).Normalized();
    Point norm(-dir.y, dir.x);
    double r = thickness / 2.0;

    Point p1 = A + norm * r;
    Point p2 = B + norm * r;
    Point p3 = B - norm * r;
    Point p4 = A - norm * r;

    if (cap == LineCap::kSquare) {
        p1 = p1 - dir * r;
        p4 = p4 - dir * r;
        p2 = p2 + dir * r;
        p3 = p3 + dir * r;
    }

    std::vector<Point2D> box = {p1.ToPoint2D(), p2.ToPoint2D(), p3.ToPoint2D(), p4.ToPoint2D()};
    FillPolygonNonZero(img, box, color);

    if (cap == LineCap::kRound) {
        FillCircle(img, A, r, color);
        FillCircle(img, B, r, color);
    }
}

static auto cross_product(const Point2D &a, const Point2D &b, const Point2D &c) -> int64_t
{
    return ((int64_t)b.x - (int64_t)a.x) * ((int64_t)c.y - (int64_t)a.y) -
           ((int64_t)b.y - (int64_t)a.y) * ((int64_t)c.x - (int64_t)a.x);
}

static auto is_on_segment(const Point2D &a, const Point2D &b, const Point2D &p) -> bool
{
    if (cross_product(a, b, p) != 0) {
        return false;
    }
    return (p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) && p.y >= std::min(a.y, b.y) &&
            p.y <= std::max(a.y, b.y));
}

static auto is_point_on_edges(const Point2D &p, const std::vector<Point2D> &poly) -> bool
{
    auto n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        if (is_on_segment(poly[i], poly[(i + 1) % n], p)) {
            return true;
        }
    }
    return false;
}

static auto segments_intersect_exact(const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d) -> bool
{
    auto cp1 = cross_product(a, b, c);
    auto cp2 = cross_product(a, b, d);
    auto cp3 = cross_product(c, d, a);
    auto cp4 = cross_product(c, d, b);

    if (((cp1 > 0 && cp2 < 0) || (cp1 < 0 && cp2 > 0)) && ((cp3 > 0 && cp4 < 0) || (cp3 < 0 && cp4 > 0))) {
        return true;
    }

    if (cp1 == 0 && is_on_segment(a, b, c)) {
        return true;
    }
    if (cp2 == 0 && is_on_segment(a, b, d)) {
        return true;
    }
    if (cp3 == 0 && is_on_segment(c, d, a)) {
        return true;
    }
    if (cp4 == 0 && is_on_segment(c, d, b)) {
        return true;
    }
    return false;
}

auto IsSimplePolygon(const std::vector<Point2D> &v) -> bool
{
    size_t n = v.size();
    if (n < 4) {
        return true;
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) {
                continue;
            }

            if (segments_intersect_exact(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n])) {
                return false;
            }
        }
    }

    return true;
}

auto IsConvexPolygon(const std::vector<Point2D> &v) -> bool
{
    size_t n = v.size();
    if (n < 3) {
        return false;
    }
    if (!IsSimplePolygon(v)) {
        return false;
    }

    auto has_pos = false;
    auto has_neg = false;

    for (size_t i = 0; i < n; ++i) {
        auto cp = cross_product(v[i], v[(i + 1) % n], v[(i + 2) % n]);
        if (cp > 0) {
            has_pos = true;
        }
        if (cp < 0) {
            has_neg = true;
        }
        if (has_pos && has_neg) {
            return false;
        }
    }
    return true;
}

auto FillPolygonEvenOdd(Image &img, const std::vector<Point2D> &vertices, Pixel color) -> void
{
    if (vertices.empty()) {
        return;
    }

    if (vertices.size() == 1) {
        if (vertices[0].x < img.WidthRgb() && vertices[0].y < img.HeightRgb()) {
            img.GetRGBPixel(vertices[0].x, vertices[0].y) = color;
        }
        return;
    }

    if (vertices.size() == 2) {
        std::ignore = DrawLine(img, vertices[0], vertices[1], color);
        return;
    }

    int min_y = img.HeightRgb(), max_y = 0;
    for (const auto &p : vertices) {
        min_y = std::min(min_y, (int)p.y);
        max_y = std::max(max_y, (int)p.y);
    }

    min_y = std::max(0, min_y);
    max_y = std::min((int)img.HeightRgb() - 1, max_y);

    for (int y = min_y; y <= max_y; ++y) {
        std::vector<int> nodes;
        auto n = vertices.size();

        for (size_t i = 0; i < n; ++i) {
            auto p1 = vertices[i];
            auto p2 = vertices[(i + 1) % n];

            if (p1.y == p2.y) {
                continue;
            }

            if ((p1.y < y && p2.y >= y) || (p2.y < y && p1.y >= y)) {
                auto num     = (double)y - (double)p1.y;
                auto den     = (double)p2.y - (double)p1.y;
                auto delta_x = (double)p2.x - (double)p1.x;

                auto x = (double)p1.x + (num / den) * delta_x;
                nodes.push_back((int)std::round(x));
            }
        }
        std::sort(nodes.begin(), nodes.end());

        for (size_t i = 0; i + 1 < nodes.size(); i += 2) {
            auto x_start = std::max(0, nodes[i]);
            auto x_end   = std::min((int)img.WidthRgb() - 1, nodes[i + 1]);

            if (x_start > x_end) {
                continue;
            }

            for (auto x = x_start; x <= x_end; ++x) {
                img.GetRGBPixel(x, y) = color;
            }
        }
    }
    DrawPolygonEdges(img, vertices, color);
}

auto FillPolygonNonZero(Image &img, const std::vector<Point2D> &poly, Pixel color) -> void
{
    if (poly.empty()) {
        return;
    }

    auto minx = poly[0].x, maxx = poly[0].x;
    auto miny = poly[0].y, maxy = poly[0].y;
    for (const auto &p : poly) {
        minx = std::min(minx, p.x);
        maxx = std::max(maxx, p.x);
        miny = std::min(miny, p.y);
        maxy = std::max(maxy, p.y);
    }
    minx = std::max(minx, size_t(0));
    miny = std::max(miny, size_t(0));
    maxx = std::min(maxx, img.WidthRgb() - 1);
    maxy = std::min(maxy, img.HeightRgb() - 1);

    for (auto y = miny; y <= maxy; y++) {
        for (auto x = minx; x <= maxx; x++) {
            Point2D p{x, y};

            if (is_point_on_edges(p, poly)) {
                img.GetRGBPixel(x, y) = color;
                continue;
            }

            auto wn = 0;
            for (size_t i = 0; i < poly.size(); ++i) {
                auto v1 = poly[i];
                auto v2 = poly[(i + 1) % poly.size()];

                if (v1.y <= p.y) {
                    if (v2.y > p.y) {
                        // an upward crossing
                        if (cross_product(v1, v2, p) > 0) { // P left of edge
                            ++wn;
                        }
                    }
                } else {
                    if (v2.y <= p.y) {
                        // a downward crossing
                        if (cross_product(v1, v2, p) < 0) { // P right of edge
                            --wn;
                        }
                    }
                }
            }

            if (wn != 0) {
                img.GetRGBPixel(x, y) = color;
            }
        }
    }
}

static double getPolygonArea(const std::vector<Point> &polygon)
{
    double area = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Point &p1 = polygon[i];
        const Point &p2 = polygon[(i + 1) % polygon.size()];
        area += (p2.x - p1.x) * (p2.y + p1.y);
    }
    return area;
}

static bool isOnSegment(Point a, Point b, Point p)
{
    double cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    if (std::abs(cross) > Point::EPS) {
        return false;
    }

    return p.x >= std::min(a.x, b.x) - Point::EPS && p.x <= std::max(a.x, b.x) + Point::EPS &&
           p.y >= std::min(a.y, b.y) - Point::EPS && p.y <= std::max(a.y, b.y) + Point::EPS;
}

static bool intersectSegments(Point &p1, Point &p2, Point v1, Point v2)
{
    Point dP = p2 - p1;
    Point dV = v2 - v1;

    double det = dP.x * dV.y - dP.y * dV.x;

    if (std::abs(det) > 1e-9) {
        double t = ((v1.x - p1.x) * dV.y - (v1.y - p1.y) * dV.x) / det;
        double u = ((v1.x - p1.x) * dP.y - (v1.y - p1.y) * dP.x) / det;

        if (t >= -1e-9 && t <= 1.000000001 && u >= -1e-9 && u <= 1.000000001) {
            p1 = p1 + dP * t;
            p2 = p1;
            return true;
        }
    } else {
        if (std::abs((v1.x - p1.x) * dP.y - (v1.y - p1.y) * dP.x) > 1e-9) {
            return false;
        }

        auto clip1D = [](double minP, double maxP, double minV, double maxV, double &resMin, double &resMax) {
            resMin = std::max(minP, minV);
            resMax = std::min(maxP, maxV);
            return resMin <= resMax + 1e-9;
        };

        double interMinX, interMaxX, interMinY, interMaxY;
        bool intersect = clip1D(std::min(p1.x, p2.x),
                                std::max(p1.x, p2.x),
                                std::min(v1.x, v2.x),
                                std::max(v1.x, v2.x),
                                interMinX,
                                interMaxX) &&
                         clip1D(std::min(p1.y, p2.y),
                                std::max(p1.y, p2.y),
                                std::min(v1.y, v2.y),
                                std::max(v1.y, v2.y),
                                interMinY,
                                interMaxY);

        if (intersect) {
            Point start = p1, end = p2;
            double tStart = (std::abs(dP.x) > std::abs(dP.y)) ? (interMinX - start.x) / dP.x
                                                              : (interMinY - start.y) / dP.y;
            double tEnd   = (std::abs(dP.x) > std::abs(dP.y)) ? (interMaxX - start.x) / dP.x
                                                              : (interMaxY - start.y) / dP.y;

            p1 = start + dP * std::clamp(std::min(tStart, tEnd), 0.0, 1.0);
            p2 = start + dP * std::clamp(std::max(tStart, tEnd), 0.0, 1.0);
            return true;
        }
    }
    return false;
}

auto CyrusBeckClipSegment(Point &p0, Point &p1, const std::vector<Point> &clipPoly) -> bool
{
    int n = static_cast<int>(clipPoly.size());
    if (n == 0) {
        return false;
    }

    double area = getPolygonArea(clipPoly);

    if (std::abs(area) < 1e-9) {
        Point polyMin = clipPoly[0];
        Point polyMax = clipPoly[0];

        for (const auto &pt : clipPoly) {
            if (pt.x < polyMin.x || (pt.x == polyMin.x && pt.y < polyMin.y)) {
                polyMin = pt;
            }
            if (pt.x > polyMax.x || (pt.x == polyMax.x && pt.y > polyMax.y)) {
                polyMax = pt;
            }
        }

        if (polyMin.Equal(polyMax)) {
            if (isOnSegment(p0, p1, polyMin)) {
                p0 = p1 = polyMin;
                return true;
            }
            return false;
        }

        return intersectSegments(p0, p1, polyMin, polyMax);
    }

    Point center{0, 0};
    for (const auto &p : clipPoly) {
        center = center + p;
    }
    center = center / (double)clipPoly.size();

    auto t_enter = 0.0;
    auto t_leave = 1.0;
    auto D       = p1 - p0;

    for (size_t i = 0; i < clipPoly.size(); ++i) {
        auto p_curr = clipPoly[i];
        auto p_next = clipPoly[(i + 1) % clipPoly.size()];

        auto edge = p_next - p_curr;
        Point normal{-edge.y, edge.x};

        auto to_edge = (p_curr + p_next) * 0.5 - center;
        if (normal.Dot(to_edge) < 0) {
            normal = normal * -1.0;
        }

        auto num = normal.Dot(p_curr - p0);
        auto den = normal.Dot(D);

        if (std::abs(den) < Point::EPS) {
            if (num < 0) {
                return false;
            }
        } else {
            double t = num / den;
            if (den > 0) {
                t_leave = std::min(t_leave, t);
            } else {
                t_enter = std::max(t_enter, t);
            }
        }
    }

    if (t_enter > t_leave) {
        return false;
    }

    auto saved_p0 = p0;
    p0            = saved_p0 + D * t_enter;
    p1            = saved_p0 + D * t_leave;
    return true;
}

static auto get_signed_distance(const Point &p, const Point &a, const Point &b) -> double
{
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
}

static auto intersect_lines_stable(const Point &s, const Point &e, const Point &a, const Point &b) -> Point
{
    auto d1 = get_signed_distance(s, a, b);
    auto d2 = get_signed_distance(e, a, b);

    if (std::abs(d1 - d2) < 1e-12) {
        return s;
    }

    auto t = d1 / (d1 - d2);
    return s + (e - s) * t;
}

auto ClipPolygonSutherlandHodgman(const std::vector<Point> &subjectPoly, const std::vector<Point> &clipPoly)
    -> std::vector<Point>
{
    double area = getPolygonArea(clipPoly);

    if (clipPoly.size() < 3) {
        auto pMin = clipPoly[0], pMax = clipPoly[0];
        for (const auto &p : clipPoly) {
            if (p.x < pMin.x || (p.x == pMin.x && p.y < pMin.y)) {
                pMin = p;
            }
            if (p.x > pMax.x || (p.x == pMax.x && p.y > pMax.y)) {
                pMax = p;
            }
        }

        std::vector<Point> result;
        if (pMin.Equal(pMax)) { // point
            for (size_t i = 0; i < subjectPoly.size(); ++i) {
                if (isOnSegment(subjectPoly[i], subjectPoly[(i + 1) % subjectPoly.size()], pMin)) {
                    result.push_back(pMin);
                    break;
                }
            }
        } else { // line
            for (auto i = 0u; i < subjectPoly.size(); ++i) {
                auto v1 = subjectPoly[i];
                auto v2 = subjectPoly[(i + 1) % subjectPoly.size()];

                auto p1 = pMin, p2 = pMax;
                if (intersectSegments(p1, p2, v1, v2)) {
                    result.push_back(p1);
                    if (!p1.Equal(p2)) {
                        result.push_back(p2);
                    }
                }
            }

            if (result.size() > 1) {
                std::sort(result.begin(), result.end(), [](const Point &a, const Point &b) {
                    return (a.x != b.x) ? a.x < b.x : a.y < b.y;
                });
                result.erase(std::unique(result.begin(),
                                         result.end(),
                                         [](const Point &a, const Point &b) { return a.Equal(b); }),
                             result.end());
            }
        }
        return result;
    }

    if (std::abs(area) < 1e-9) {
        auto pMin = clipPoly[0], pMax = clipPoly[0];
        for (const auto &p : clipPoly) {
            if (p.x < pMin.x || (p.x == pMin.x && p.y < pMin.y)) {
                pMin = p;
            }
            if (p.x > pMax.x || (p.x == pMax.x && p.y > pMax.y)) {
                pMax = p;
            }
        }

        std::vector<Point> result;
        if (pMin.Equal(pMax)) {
            return result;
        }

        for (size_t i = 0; i < subjectPoly.size(); ++i) {
            Point a     = subjectPoly[i];
            Point b     = subjectPoly[(i + 1) % subjectPoly.size()];
            Point start = a, end = b;
            if (intersectSegments(start, end, pMin, pMax)) {
                result.push_back(start);
                if (!start.Equal(end)) {
                    result.push_back(end);
                }
            }
        }
        result.erase(
            std::unique(result.begin(), result.end(), [](const Point &a, const Point &b) { return a.Equal(b); }),
            result.end());
        return result;
    }

    auto is_ccw = (area > 0);

    auto output = subjectPoly;

    for (size_t i = 0; i < clipPoly.size(); ++i) {
        auto a = clipPoly[i];
        auto b = clipPoly[(i + 1) % clipPoly.size()];

        auto input = output;
        output.clear();
        if (input.empty()) {
            break;
        }

        auto s = input.back();
        for (const auto &e : input) {
            auto dist_e = get_signed_distance(e, a, b);
            auto dist_s = get_signed_distance(s, a, b);

            auto is_inside = [&](double dist) { return is_ccw ? (dist >= -1e-9) : (dist <= 1e-9); };

            if (is_inside(dist_e)) {
                if (!is_inside(dist_s)) {
                    output.push_back(intersect_lines_stable(s, e, a, b));
                }
                output.push_back(e);
            } else if (is_inside(dist_s)) {
                output.push_back(intersect_lines_stable(s, e, a, b));
            }
            s = e;
        }
    }
    return output;
}

auto BezierCubicCurve(Point p0, Point p1, Point p2, Point p3) -> std::vector<Point>
{
    auto len = (p1 - p0).Len() + (p2 - p1).Len() + (p3 - p2).Len();

    auto steps = std::max(10, static_cast<int>(len * 2.0));
    if (steps > 5000) {
        steps = 5000;
    }

    std::vector<Point> curve;
    curve.reserve(steps + 1);

    for (auto i = 0; i <= steps; i++) {
        auto t   = double(i) / steps;
        auto u   = 1.0 - t;
        auto tt  = t * t;
        auto uu  = u * u;
        auto uuu = uu * u;
        auto ttt = tt * t;

        Point p = p0 * uuu + p1 * (3 * uu * t) + p2 * (3 * u * tt) + p3 * ttt;
        curve.push_back(p);
    }
    return curve;
}

static auto rotate_point(Point3 p, Point3 axis, double angle) -> Point3
{
    axis           = axis.Normalized();
    auto cos_theta = std::cos(angle);
    auto sin_theta = std::sin(angle);

    Point3 cross = axis.Cross(p);
    auto dot     = axis.Dot(p);

    return p * cos_theta + cross * sin_theta + axis * dot * (1.0 - cos_theta);
}

static auto project_point(Point3 p, ProjectionType type, double k, double center_x, double center_y) -> Point
{
    double x_proj, y_proj;

    if (type == ProjectionType::kParallel) {
        x_proj = p.x;
        y_proj = p.y;
    } else {
        double dist = k - p.z;
        if (dist < 1.0) {
            dist = 1.0;
        }
        double factor = k / dist;
        x_proj        = p.x * factor;
        y_proj        = p.y * factor;
    }

    return Point(center_x + x_proj, center_y - y_proj);
}

auto RenderParallelepiped(Image &img,
                          Point3 obj_center,
                          Point3 size,
                          Point3 rotationAxis,
                          double angle,
                          ProjectionType type,
                          double k) -> void
{
    auto w = size.x / 2.0;
    auto h = size.y / 2.0;
    auto d = size.z / 2.0;

    std::vector<Point3> local_verts = {
        {-w, -h,  d},
        { w, -h,  d},
        { w,  h,  d},
        {-w,  h,  d}, // z > 0
        {-w, -h, -d},
        { w, -h, -d},
        { w,  h, -d},
        {-w,  h, -d}  // z < 0
    };

    std::vector<std::vector<int> > faces = {
        {0, 1, 2, 3}, // front
        {1, 5, 6, 2}, // right
        {5, 4, 7, 6}, // back
        {4, 0, 3, 7}, // left
        {3, 2, 6, 7}, // top
        {4, 5, 1, 0}  // bottom
    };

    Matrix4x4 rot = Matrix4x4::Rotation(rotationAxis, angle);

    std::vector<Point3> world_verts;
    for (const auto &v : local_verts) {
        Point3 worldP = rot.Transform(v) + obj_center;
        world_verts.push_back(worldP);
    }

    auto screen_cx = img.WidthRgb() / 2.0;
    auto screen_cy = img.HeightRgb() / 2.0;

    auto cameraPos = (type == ProjectionType::kPerspective) ? Point3(0, 0, k) : Point3(0, 0, 1000000);

    for (const auto &face_indices : faces) {
        const auto &p0 = world_verts[face_indices[0]];
        const auto &p1 = world_verts[face_indices[1]];
        const auto &p2 = world_verts[face_indices[2]];

        auto v1     = p1 - p0;
        auto v2     = p2 - p0;
        auto normal = v1.Cross(v2).Normalized();

        auto face_center = (p0 + p2) * 0.5;
        Point3 view_vec;

        if (type == ProjectionType::kParallel) {
            view_vec = Point3(0, 0, 1);
        } else {
            view_vec = (cameraPos - face_center).Normalized();
        }

        if (normal.Dot(view_vec) > 0) {
            std::vector<Point2D> poly2d;
            for (int idx : face_indices) {
                auto proj = project_point(world_verts[idx], type, k, screen_cx, screen_cy);
                poly2d.push_back(proj.ToPoint2D());
            }

            DrawPolygonEdges(img, poly2d, {200, 255, 200});
        }
    }
}

auto ExtractExteriorContour(const std::vector<Point2D> &poly, size_t width, size_t height) -> std::vector<Point2D>
{
    Image maskRGB(width, height, Image::Type::kRGB);

    Pixel black{0, 0, 0};
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            maskRGB.GetRGBPixel(x, y) = black;
        }
    }

    FillPolygonNonZero(maskRGB, poly, {255, 255, 255});

    auto start_x = -1, start_y = -1;
    auto found = false;

    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            if (maskRGB.GetRGBPixel(x, y).r == 255) {
                start_x = static_cast<int>(x);
                start_y = static_cast<int>(y);
                found   = true;
                break;
            }
        }
        if (found) {
            break;
        }
    }

    if (!found) {
        return {};
    }

    std::vector<Point2D> contour;
    contour.push_back({(size_t)start_x, (size_t)start_y});

    int dx[] = {0, 1, 1, 1, 0, -1, -1, -1};
    int dy[] = {-1, -1, 0, 1, 1, 1, 0, -1};

    auto cx        = start_x;
    auto cy        = start_y;
    auto backtrack = 6;

    auto max_iters = static_cast<int>(width * height * 4);
    auto iter      = 0;

    while (iter++ < max_iters) {
        auto found_neighbor = -1;

        for (auto i = 0; i < 8; ++i) {
            auto idx = (backtrack + 1 + i) % 8;
            auto nx  = cx + dx[idx];
            auto ny  = cy + dy[idx];

            if (nx >= 0 && nx < static_cast<int>(width) && ny >= 0 && ny < static_cast<int>(height)) {
                if (maskRGB.GetRGBPixel(nx, ny).r == 255) {
                    found_neighbor = idx;
                    cx             = nx;
                    cy             = ny;
                    break;
                }
            }
        }

        if (found_neighbor == -1) {
            break;
        }
        if (cx == start_x && cy == start_y) {
            break;
        }

        contour.push_back({(size_t)cx, (size_t)cy});
        backtrack = (found_neighbor + 4) % 8;
    }

    return contour;
}

static std::vector<Point> CreateSingleBezierArc(Point center, double radius, double start_rad, double end_rad)
{
    auto theta = end_rad - start_rad;
    auto k     = 4.0 / 3.0 * std::tan(theta / 4.0);

    auto p0 = center + Point(std::cos(start_rad), std::sin(start_rad)) * radius;
    auto p3 = center + Point(std::cos(end_rad), std::sin(end_rad)) * radius;

    Point t0(-std::sin(start_rad), std::cos(start_rad));
    Point t3(-std::sin(end_rad), std::cos(end_rad));

    auto p1 = p0 + t0 * (k * radius);
    auto p2 = p3 - t3 * (k * radius);

    return BezierCubicCurve(p0, p1, p2, p3);
}

auto MakeCircleArc(Point center, double radius, double angle_start_deg, double angle_end_deg) -> std::vector<Point>
{
    std::vector<Point> result;

    auto to_rad = M_PI / 180.0;
    auto start  = angle_start_deg * to_rad;
    auto end    = angle_end_deg * to_rad;

    if (end < start) {
        end += 2.0 * M_PI;
    }

    auto total_sweep    = end - start;
    const auto max_step = M_PI / 2.0;

    auto segments = static_cast<int>(std::ceil(total_sweep / max_step));
    auto step     = total_sweep / segments;

    for (auto i = 0; i < segments; ++i) {
        auto s           = start + i * step;
        auto e           = start + (i + 1) * step;
        auto arc_segment = CreateSingleBezierArc(center, radius, s, e);

        if (i > 0 && !result.empty()) {
            result.insert(result.end(), arc_segment.begin() + 1, arc_segment.end());
        } else {
            result.insert(result.end(), arc_segment.begin(), arc_segment.end());
        }
    }
    return result;
}

auto ColorQuantizationKMeans(Image &img, int k) -> std::expected<void, boost::system::error_code>
{
    struct ClusterData {
        double sum = 0.0;
        int count  = 0;
    };

    if (k <= 0) {
        return {};
    }

    const int rows        = img.HeightGray();
    const int cols        = img.WidthGray();
    const auto num_pixels = rows * cols;

    if (num_pixels == 0) {
        return {};
    }

    const auto safe_k = (static_cast<size_t>(k) > static_cast<size_t>(num_pixels)) ? num_pixels : k;

    auto max_attempt = 100;
    std::vector<double> centroids(safe_k);
    std::set<uint8_t> used_colors;

    std::mt19937 gen(42);
    std::uniform_int_distribution<int> dist_x(0, cols - 1);
    std::uniform_int_distribution<int> dist_y(0, rows - 1);

    for (auto i = 0; i < safe_k; ++i) {
        auto color = img(dist_x(gen), dist_y(gen));
        if (!used_colors.contains(color)) {
            centroids[i] = static_cast<double>(color);
            used_colors.insert(color);
        } else {
            --i;
            --max_attempt;
            if (max_attempt <= 0) {
                return std::unexpected{boost::system::error_code{}};
            }
        }
    }

    const auto max_iters = 200;
    std::vector<int> labels(num_pixels, -1);

    auto get_dist_sq = [](double val, double center) { return (val - center) * (val - center); };

    for (auto iter = 0; iter < max_iters; ++iter) {
        std::vector<ClusterData> new_centroids(safe_k);

        auto changed = false;

        for (auto y = 0; y < rows; ++y) {
            for (auto x = 0; x < cols; ++x) {
                auto val = img(x, y);

                auto min_d  = std::numeric_limits<double>::max();
                auto best_c = 0;

                for (auto j = 0; j < safe_k; ++j) {
                    auto d = get_dist_sq(val, centroids[j]);
                    if (d < min_d) {
                        min_d  = d;
                        best_c = j;
                    }
                }

                auto linear_idx = y * cols + x;
                if (labels[linear_idx] != best_c) {
                    changed            = true;
                    labels[linear_idx] = best_c;
                }

                new_centroids[best_c].sum += val;
                new_centroids[best_c].count++;
            }
        }

        for (auto j = 0; j < safe_k; ++j) {
            if (new_centroids[j].count > 0) {
                centroids[j] = new_centroids[j].sum / new_centroids[j].count;
            } else {
                centroids[j] = static_cast<double>(img(dist_x(gen), dist_y(gen)));
            }
        }

        if (!changed) {
            break;
        }
    }

    for (auto y = 0; y < rows; ++y) {
        for (auto x = 0; x < cols; ++x) {
            auto linear_idx = y * cols + x;
            auto c          = labels[linear_idx];

            img(x, y) = static_cast<Image::gray_type>(std::clamp(centroids[c], 0.0, 255.0));
        }
    }

    return {};
}
} // namespace improcessing
