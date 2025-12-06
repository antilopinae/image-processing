#include <cmath>
#include <cstdlib>
#include <file_closer.hpp>
#include <image.hpp>
#include <improcessing.hpp>
#include <png_guard.hpp>
#include <string>
#include <vector>

namespace improcessing {
    auto ReadImage(const std::string &filename) -> std::expected<Image, boost::system::error_code> {
        std::unique_ptr<FILE, details::FileCloser> file(std::fopen(filename.c_str(), "rb"));
        if (!file) {
            return std::unexpected{
                boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
            };
        }

        PngGuard png(false);

        png_init_io(png.png_ptr, file.get());
        png_read_info(png.png_ptr, png.info_ptr);

        png_uint_32 WidthGray = 0;
        png_uint_32 HeightGray = 0;
        int bit_depth = 0;
        int color_type = 0;

        png_get_IHDR(png.png_ptr, png.info_ptr, &WidthGray, &HeightGray, &bit_depth, &color_type, NULL, NULL, NULL);

        if (color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_palette_to_rgb(png.png_ptr);

        if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
            png_set_rgb_to_gray_fixed(png.png_ptr, 1, -1, -1);

        if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_strip_alpha(png.png_ptr);

        if (bit_depth == 16)
            png_set_strip_16(png.png_ptr);

        png_read_update_info(png.png_ptr, png.info_ptr);

        Image image(WidthGray, HeightGray);

        auto row_ptrs = image.GetGrayPtrs();

        png_read_image(png.png_ptr, row_ptrs.data());
        png_read_end(png.png_ptr, png.info_ptr);

        return std::move(image);
    }

    auto SaveImage(const std::string &filename, const Image &image,
                   Image::Type type) -> std::expected<void, boost::system::error_code> {
        std::unique_ptr<FILE, details::FileCloser> file(std::fopen(filename.c_str(), "wb"));
        if (!file) {
            return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::io_error)};
        }

        PngGuard png(true);

        png_init_io(png.png_ptr, file.get());

        switch (type) {
            case Image::Type::kGray: {
                png_set_IHDR(
                    png.png_ptr,
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
            }
            break;
            case Image::Type::kRGB: {
                png_set_IHDR(
                    png.png_ptr,
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
            }
            break;
        }

        png_write_end(png.png_ptr, png.info_ptr);

        return {};
    }

    auto MakeCircularGrayscale(Image &image, double radius_fraction) -> std::expected<void, boost::system::error_code> {
        const auto cx = (image.WidthGray() - 1) / 2.0;
        const auto cy = (image.HeightGray() - 1) / 2.0;
        const auto r = std::min(image.WidthGray(), image.HeightGray()) * radius_fraction;
        const auto edge_WidthGray = std::max(1.0, r * 0.02);

        for (int y = 0; y < image.HeightGray(); ++y) {
            for (int x = 0; x < image.WidthGray(); ++x) {
                const auto dx = x - cx;
                const auto dy = y - cy;
                const auto d = std::sqrt(dx * dx + dy * dy);

                uint8_t value = 0;

                if (d <= r - edge_WidthGray) {
                    const auto t = d / (r - edge_WidthGray);
                    auto v = 255.0 * (1.0 - 0.7 * t);

                    v = std::clamp(v, 0.0, 255.0);
                    value = static_cast<uint8_t>(std::round(v));
                } else if (d <= r + edge_WidthGray) {
                    const auto t = (d - (r - edge_WidthGray)) / (2.0 * edge_WidthGray);
                    auto v = 255.0 * (1.0 - t);

                    v = std::clamp(v, 0.0, 255.0);
                    value = static_cast<uint8_t>(std::round(v));
                }

                image.DataGray()[y * image.WidthGray() + x] = value;
            }
        }

        return {};
    }

    auto Blend(const Image &A, const Image &B, const Image &Alpha) -> std::expected<Image, boost::system::error_code> {
        if (A.WidthGray() > B.WidthGray() || A.HeightGray() > B.HeightGray() || A.WidthGray() > Alpha.WidthGray() || A.
            HeightGray() > Alpha.
            HeightGray()) {
            return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::invalid_argument)};
        }

        Image out;
        out.ResizeGray(A.WidthGray(), A.HeightGray());

        const auto total = A.WidthGray() * A.HeightGray();

        for (size_t i = 0; i < total; ++i) {
            const auto a = A.DataGray()[i];
            const auto b = B.DataGray()[i];
            const auto alpha = Alpha.DataGray()[i];

            const auto blended = static_cast<uint64_t>(alpha) * b + static_cast<uint64_t>(255 - alpha) * a;

            out.DataGray()[i] = static_cast<uint8_t>(blended / 255);
        }

        return std::move(out);
    }

    auto FloydSteinbergDither(Image &img, uint8_t n_levels) -> std::expected<void, boost::system::error_code> {
        const int rows = img.HeightGray();
        const int cols = img.WidthGray();
        const int levels = 1 << n_levels;

        struct {
            Image::gray_type newval;
        } lut[256];

        for (int i = 0; i < 256; ++i) {
            int nearest = (i * (levels - 1) + 127) / 255;
            int quantized = (nearest * 255 + (levels - 1) / 2) / (levels - 1);
            lut[i] = {static_cast<Image::gray_type>(quantized)};
        }

        int err;
        auto apply_error = [&](int x, int y, int num, int den) {
            if (y >= 0 && y < rows && x >= 0 && x < cols) {
                int val = img(x, y) + (err * num) / den;
                img(x, y) = std::clamp(val, 0, 255);
            }
        };

        for (int y = 0; y < rows; ++y) {
            const bool reverse = y & 1;
            const auto start = reverse ? cols - 1 : 0;
            const auto end = reverse ? -1 : cols;
            const auto step = reverse ? -1 : 1;

            for (int x = start; x != end; x += step) {
                const auto old_val = img(x, y);
                const auto entry = lut[old_val];
                const auto new_val = entry.newval;
                img(x, y) = new_val;

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

    auto DrawLine(Image &img, Point2D start, Point2D end,
                  Pixel color) -> std::expected<void, boost::system::error_code> {
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

    auto DrawPolygonEdges(Image &img, const std::vector<Point2D> &poly,
                          Pixel color) -> std::expected<void, boost::system::error_code> {
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

    __always_inline static auto orient(const Point2D &a, const Point2D &b, const Point2D &c) -> int64_t {
        int64_t abx = b.x - a.x;
        int64_t aby = b.y - a.y;
        int64_t acx = c.x - a.x;
        int64_t acy = c.y - a.y;
        return abx * acy - aby * acx;
    }

    __always_inline static auto onSegment(const Point2D &a, const Point2D &b, const Point2D &p) -> bool {
        if (__builtin_expect(orient(a, b, p) != 0, 1)) return false;

        return (p.x >= (a.x < b.x ? a.x : b.x) &&
                p.x <= (a.x > b.x ? a.x : b.x) &&
                p.y >= (a.y < b.y ? a.y : b.y) &&
                p.y <= (a.y > b.y ? a.y : b.y));
    }

    __always_inline static auto segmentsIntersect(
        const Point2D &p1, const Point2D &p2,
        const Point2D &q1, const Point2D &q2) -> bool {
        // AABB
        if (std::max(p1.x, p2.x) < std::min(q1.x, q2.x)) return false;
        if (std::max(q1.x, q2.x) < std::min(p1.x, p2.x)) return false;
        if (std::max(p1.y, p2.y) < std::min(q1.y, q2.y)) return false;
        if (std::max(q1.y, q2.y) < std::min(p1.y, p2.y)) return false;

        int64_t o1 = orient(p1, p2, q1);
        int64_t o2 = orient(p1, p2, q2);
        if (__builtin_expect((int64_t) o1 * o2 < 0, 0)) {
            int64_t o3 = orient(q1, q2, p1);
            int64_t o4 = orient(q1, q2, p2);
            return (o3 * o4 < 0);
        }

        if (o1 == 0 && onSegment(p1, p2, q1)) return true;
        if (o2 == 0 && onSegment(p1, p2, q2)) return true;

        int64_t o3 = orient(q1, q2, p1);
        int64_t o4 = orient(q1, q2, p2);
        if (o3 == 0 && onSegment(q1, q2, p1)) return true;
        if (o4 == 0 && onSegment(q1, q2, p2)) return true;

        return false;
    }

    auto IsSimplePolygon(const std::vector<Point2D> &v) -> bool {
        auto n = v.size();
        if (n < 4) return true;

        for (auto i = 0ul; i < n; i++) {
            auto i2 = (i + 1 < n ? i + 1 : 0ul);
            for (auto j = i + 1; j < n; j++) {
                auto j2 = (j + 1 < n ? j + 1 : 0ul);

                if (i == j) continue;
                if (i2 == j) continue;
                if (j2 == i) continue;

                if (segmentsIntersect(v[i], v[i2], v[j], v[j2]))
                    return false;
            }
        }
        return true;
    }

    auto IsConvexPolygon(const std::vector<Point2D> &v) -> bool {
        auto n = v.size();
        if (n < 4) return true;

        int s0 = 0;

        for (auto i = 0ul; i < n; i++) {
            auto i1 = (i + 1 < n ? i + 1 : 0ul);
            auto i2 = (i1 + 1 < n ? i1 + 1 : 0ul);

            int64_t cr = orient(v[i], v[i1], v[i2]);
            if (cr == 0) continue;

            int s = (cr > 0 ? 1 : -1);

            if (s0 == 0) s0 = s;
            else if (s0 != s) return false;
        }
        return true;
    }

    __always_inline static auto pointOnEdge(const Point2D &p,
                                            const std::vector<Point2D> &v) -> bool {
        const auto n = v.size();
        for (auto i = 0ul; i < n; i++) {
            auto j = (i == n - 1 ? 0ul : i + 1);

            if (onSegment(v[i], v[j], p)) {
                return true;
            }
        }
        return false;
    }

    __always_inline static auto pointInPolygonNonZero(const Point2D &P, const std::vector<Point2D> &poly) {
        if (pointOnEdge(P, poly)) return true;

        int winding = 0;
        int n = (int) poly.size();
        for (int i = 0; i < n; ++i) {
            const Point2D &a = poly[i];
            const Point2D &b = poly[(i + 1 == n) ? 0 : i + 1];

            if (a.y == b.y) continue;

            if (a.y <= P.y && b.y > P.y) {
                if (orient(a, b, P) > 0)
                    ++winding;
            } else if (a.y > P.y && b.y <= P.y) {
                if (orient(a, b, P) < 0)
                    --winding;
            }
        }
        return winding != 0;
    }

    auto FillPolygonNonZero(Image &img,
                            const std::vector<Point2D> &poly,
                            Pixel color) -> void {
        auto minx = poly[0].x, maxx = poly[0].x;
        auto miny = poly[0].y, maxy = poly[0].y;

        for (auto &p: poly) {
            minx = std::min(minx, p.x);
            maxx = std::max(maxx, p.x);
            miny = std::min(miny, p.y);
            maxy = std::max(maxy, p.y);
        }

        for (auto y = miny; y <= maxy; y++)
            for (auto x = minx; x <= maxx; x++)
                if (pointInPolygonNonZero({x, y}, poly))
                    img.GetRGBPixel(x, y) = color;
    }

    __always_inline static auto findIntersectionsY_EO(const std::vector<Point2D> &vertices,
                                                      size_t y) -> std::vector<int> {
        using namespace std;

        vector<int> xInts;
        size_t n = vertices.size();

        for (size_t i = 0; i < n; ++i) {
            auto p1 = vertices[i];
            auto p2 = vertices[(i + 1) % n];

            if (p1.y == p2.y) continue;

            if (y >= min(p1.y, p2.y) && y < max(p1.y, p2.y)) {
                int dy = p2.y - p1.y;
                int dx = p2.x - p1.x;
                int x = p1.x + (y - p1.y) * dx / dy;
                xInts.push_back(x);
            }
        }

        sort(xInts.begin(), xInts.end());
        return xInts;
    }

    auto FillPolygonEvenOdd(Image &img, const std::vector<Point2D> &vertices, Pixel color) -> void {
        if (vertices.size() < 3) return;

        using namespace std;

        int minY = (int) img.HeightRgb(), maxY = 0;
        for (const auto &p: vertices) {
            minY = min(minY, (int) p.y);
            maxY = max(maxY, (int) p.y);
        }

        for (int y = max(minY, 0); y < min(maxY, (int) img.HeightRgb()); ++y) {
            vector<int> xInts = findIntersectionsY_EO(vertices, y);

            for (size_t i = 0; i + 1 < xInts.size(); i += 2) {
                int xStart = max(xInts[i], 0);
                int xEnd = min(xInts[i + 1], (int) img.WidthRgb() - 1);
                for (int x = xStart; x <= xEnd; ++x) {
                    img.GetRGBPixel(x, y) = color;
                }
            }
        }
    }

    auto CyrusBeckClipSegmentCW(Point &p0, Point &p1, const std::vector<Point> &polyCW) -> bool {
        auto n = polyCW.size();
        if (n < 3) return false;

        auto t_enter = 0.0;
        auto t_leave = 1.0;

        auto sx = p1.x - p0.x;
        auto sy = p1.y - p0.y;

        for (size_t i = 0; i < n; ++i) {
            const auto &vi = polyCW[i];
            const auto &vi1 = polyCW[(i + 1) % n];

            auto nx = vi1.y - vi.y;
            auto ny = vi.x - vi1.x;

            auto denom = nx * sx + ny * sy;
            auto num = nx * (p0.x - vi.x) + ny * (p0.y - vi.y);

            if (fabs(denom) > Point::EPS) {
                auto t = -num / denom;
                if (denom > 0.0) {
                    if (t > t_enter) {
                        t_enter = t;
                    }
                } else {
                    if (t < t_leave) {
                        t_leave = t;
                    }
                }
                if (t_enter - t_leave > Point::EPS) {
                    return false;
                }
            } else {
                int cls = ClassifyPointToEdge(vi, vi1, p0);
                if (cls > 0) {
                    // LEFT
                    return false;
                }
            }
        }

        auto t0 = std::max(0.0, t_enter);
        auto t1 = std::min(1.0, t_leave);
        if (t0 > t1 + Point::EPS) {
            return false;
        }

        auto d = p1 - p0;
        auto c0 = p0 + d * t0;
        auto c1 = p0 + d * t1;

        // Update points
        p0 = c0;
        p1 = c1;

        return true;
    }

    auto BezierLine(const Point &p0, const Point &p1, double t) -> Point {
        return p0 * (1.0 - t) + p1 * t;
    }

    auto BezierQuadratic(const Point &p0, const Point &p1, const Point &p2, double t) -> Point {
        return BezierLine(
            BezierLine(p0, p1, t),
            BezierLine(p1, p2, t),
            t
        );
    }

    auto BezierCubic(const Point &p0, const Point &p1, const Point &p2, const Point &p3, double t) -> Point {
        return BezierLine(
            BezierQuadratic(p0, p1, p2, t),
            BezierQuadratic(p1, p2, p3, t),
            t
        );
    }

    auto BezierCubicCurve(Point p0, Point p1, Point p2, Point p3, int steps) -> std::vector<Point> {
        std::vector<Point> curve;
        curve.reserve(steps + 1);

        for (int i = 0; i <= steps; i++) {
            double t = double(i) / steps;
            curve.push_back(BezierCubic(p0, p1, p2, p3, t));
        }
        return curve;
    }
} // namespace improcessing
