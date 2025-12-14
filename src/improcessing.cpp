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

    static auto cross_product(const Point2D &a, const Point2D &b, const Point2D &c) -> int64_t {
        return ((int64_t) b.x - (int64_t) a.x) * ((int64_t) c.y - (int64_t) a.y) -
               ((int64_t) b.y - (int64_t) a.y) * ((int64_t) c.x - (int64_t) a.x);
    }

    static auto is_on_segment(const Point2D &a, const Point2D &b, const Point2D &p) -> bool {
        if (cross_product(a, b, p) != 0) return false;
        return (p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) &&
                p.y >= std::min(a.y, b.y) && p.y <= std::max(a.y, b.y));
    }

    static auto is_point_on_edges(const Point2D &p, const std::vector<Point2D> &poly) -> bool {
        size_t n = poly.size();
        for (size_t i = 0; i < n; ++i) {
            if (is_on_segment(poly[i], poly[(i + 1) % n], p)) return true;
        }
        return false;
    }

    static bool segments_intersect_exact(const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d) {
        auto cp1 = cross_product(a, b, c);
        auto cp2 = cross_product(a, b, d);
        auto cp3 = cross_product(c, d, a);
        auto cp4 = cross_product(c, d, b);

        if (((cp1 > 0 && cp2 < 0) || (cp1 < 0 && cp2 > 0)) &&
            ((cp3 > 0 && cp4 < 0) || (cp3 < 0 && cp4 > 0)))
            return true;

        if (cp1 == 0 && is_on_segment(a, b, c)) return true;
        if (cp2 == 0 && is_on_segment(a, b, d)) return true;
        if (cp3 == 0 && is_on_segment(c, d, a)) return true;
        if (cp4 == 0 && is_on_segment(c, d, b)) return true;
        return false;
    }

    auto IsSimplePolygon(const std::vector<Point2D> &v) -> bool {
        size_t n = v.size();
        if (n < 4) return true;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1) continue;

                if (segments_intersect_exact(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n])) {
                    return false;
                }
            }
        }
        return true;
    }

    auto IsConvexPolygon(const std::vector<Point2D> &v) -> bool {
        size_t n = v.size();
        if (n < 3) return false;
        if (!IsSimplePolygon(v)) return false;

        bool has_pos = false;
        bool has_neg = false;

        for (size_t i = 0; i < n; ++i) {
            auto cp = cross_product(v[i], v[(i + 1) % n], v[(i + 2) % n]);
            if (cp > 0) has_pos = true;
            if (cp < 0) has_neg = true;
            if (has_pos && has_neg) return false;
        }
        return true;
    }

    auto FillPolygonEvenOdd(Image &img, const std::vector<Point2D> &vertices, Pixel color) -> void {
        if (vertices.size() < 3) return;

        int min_y = img.HeightRgb(), max_y = 0;
        for (const auto &p: vertices) {
            min_y = std::min(min_y, (int) p.y);
            max_y = std::max(max_y, (int) p.y);
        }
        min_y = std::max(0, min_y);
        max_y = std::min((int) img.HeightRgb() - 1, max_y);

        for (int y = min_y; y <= max_y; ++y) {
            std::vector<int> nodes;
            size_t n = vertices.size();
            for (size_t i = 0; i < n; ++i) {
                Point2D p1 = vertices[i];
                Point2D p2 = vertices[(i + 1) % n];

                if (p1.y == p2.y) continue;

                if ((p1.y < y && p2.y >= y) || (p2.y < y && p1.y >= y)) {
                    double num = (double) y - (double) p1.y;
                    double den = (double) p2.y - (double) p1.y;
                    double delta_x = (double) p2.x - (double) p1.x;

                    double x = (double) p1.x + (num / den) * delta_x;
                    nodes.push_back((int) std::round(x));
                }
            }
            std::sort(nodes.begin(), nodes.end());

            for (size_t i = 0; i + 1 < nodes.size(); i += 2) {
                int x_start = std::max(0, nodes[i]);
                int x_end = std::min((int) img.WidthRgb() - 1, nodes[i + 1]);

                if (x_start > x_end) continue;

                for (int x = x_start; x <= x_end; ++x) {
                    img.GetRGBPixel(x, y) = color;
                }
            }
        }
        DrawPolygonEdges(img, vertices, color);
    }

    auto FillPolygonNonZero(Image &img, const std::vector<Point2D> &poly, Pixel color) -> void {
        if (poly.empty()) return;

        auto minx = poly[0].x, maxx = poly[0].x;
        auto miny = poly[0].y, maxy = poly[0].y;
        for (const auto &p: poly) {
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

                int wn = 0;
                for (size_t i = 0; i < poly.size(); ++i) {
                    Point2D v1 = poly[i];
                    Point2D v2 = poly[(i + 1) % poly.size()];

                    if (v1.y <= p.y) {
                        if (v2.y > p.y) {
                            // an upward crossing
                            if (cross_product(v1, v2, p) > 0) // P left of edge
                                ++wn;
                        }
                    } else {
                        if (v2.y <= p.y) {
                            // a downward crossing
                            if (cross_product(v1, v2, p) < 0) // P right of edge
                                --wn;
                        }
                    }
                }

                if (wn != 0) {
                    img.GetRGBPixel(x, y) = color;
                }
            }
        }
    }

    auto CyrusBeckClipSegment(Point &p0, Point &p1, const std::vector<Point> &clipPoly) -> bool {
        if (clipPoly.size() < 3) return false;

        Point center{0, 0};
        for (const auto &p: clipPoly) center = center + p;
        center = center / (double) clipPoly.size();

        double t_enter = 0.0;
        double t_leave = 1.0;
        Point D = p1 - p0;

        for (size_t i = 0; i < clipPoly.size(); ++i) {
            Point p_curr = clipPoly[i];
            Point p_next = clipPoly[(i + 1) % clipPoly.size()];

            Point edge = p_next - p_curr;
            Point normal{-edge.y, edge.x};

            Point to_edge = (p_curr + p_next) * 0.5 - center;
            if (normal.Dot(to_edge) < 0) {
                normal = normal * -1.0;
            }

            double num = normal.Dot(p_curr - p0);
            double den = normal.Dot(D);

            if (std::abs(den) < Point::EPS) {
                if (num < 0) return false;
            } else {
                double t = num / den;
                if (den > 0) {
                    t_leave = std::min(t_leave, t);
                } else {
                    t_enter = std::max(t_enter, t);
                }
            }
        }

        if (t_enter > t_leave) return false;

        Point saved_p0 = p0;
        p0 = saved_p0 + D * t_enter;
        p1 = saved_p0 + D * t_leave;
        return true;
    }

    static double get_signed_distance(const Point &p, const Point &a, const Point &b) {
        return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    }

    static Point intersect_lines_stable(const Point &s, const Point &e, const Point &a, const Point &b) {
        double d1 = get_signed_distance(s, a, b);
        double d2 = get_signed_distance(e, a, b);

        if (std::abs(d1 - d2) < 1e-12) return s;

        double t = d1 / (d1 - d2);
        return s + (e - s) * t;
    }

    auto ClipPolygonSutherlandHodgman(const std::vector<Point> &subjectPoly,
                                      const std::vector<Point> &clipPoly) -> std::vector<Point> {
        if (clipPoly.size() < 3 || subjectPoly.empty()) return subjectPoly;

        double area = 0;
        for (size_t i = 0; i < clipPoly.size(); ++i) {
            area += clipPoly[i].x * clipPoly[(i + 1) % clipPoly.size()].y;
            area -= clipPoly[i].y * clipPoly[(i + 1) % clipPoly.size()].x;
        }
        bool is_ccw = (area > 0);

        std::vector<Point> output = subjectPoly;

        for (size_t i = 0; i < clipPoly.size(); ++i) {
            Point a = clipPoly[i];
            Point b = clipPoly[(i + 1) % clipPoly.size()];

            std::vector<Point> input = output;
            output.clear();
            if (input.empty()) break;

            Point s = input.back();
            for (const auto &e: input) {
                double dist_e = get_signed_distance(e, a, b);
                double dist_s = get_signed_distance(s, a, b);

                auto is_inside = [&](double dist) {
                    return is_ccw ? (dist >= -1e-9) : (dist <= 1e-9);
                };

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

    auto BezierCubicCurve(Point p0, Point p1, Point p2, Point p3) -> std::vector<Point> {
        double len = (p1 - p0).Len() + (p2 - p1).Len() + (p3 - p2).Len();

        int steps = std::max(10, static_cast<int>(len * 2.0));
        if (steps > 5000) steps = 5000;

        std::vector<Point> curve;
        curve.reserve(steps + 1);

        for (int i = 0; i <= steps; i++) {
            double t = double(i) / steps;
            double u = 1.0 - t;
            double tt = t * t;
            double uu = u * u;
            double uuu = uu * u;
            double ttt = tt * t;

            Point p = p0 * uuu + p1 * (3 * uu * t) + p2 * (3 * u * tt) + p3 * ttt;
            curve.push_back(p);
        }
        return curve;
    }
} // namespace improcessing
