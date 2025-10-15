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

    auto SaveImage(const std::string &filename, const Image &image) -> std::expected<void, boost::system::error_code> {
        std::unique_ptr<FILE, details::FileCloser> file(std::fopen(filename.c_str(), "wb"));
        if (!file) {
            return std::unexpected{boost::system::errc::make_error_code(boost::system::errc::io_error)};
        }

        PngGuard png(true);

        png_init_io(png.png_ptr, file.get());

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

        int err;
        auto apply_error = [&](int x, int y, int num, int den) {
            if (y >= 0 && y < rows && x >= 0 && x < cols) {
                int val = img(x, y) + (err * num) / den;
                img(x, y) = std::clamp(val, 0, 255);
            }
        };

        int levels = 1 << n_levels;

        bool reverse = false;
        int start, end, step, oldval, newval, nearest_level;

        for (int y = 0; y < rows; ++y) {
            reverse = (y % 2 != 0);
            start = reverse ? cols - 1 : 0;
            end = reverse ? -1 : cols;
            step = reverse ? -1 : 1;

            for (int x = start; x != end; x += step) {
                oldval = img(x, y);
                nearest_level = (oldval * (levels - 1) + 127) / 255;
                newval = (nearest_level * 255 + (levels - 1) / 2) / (levels - 1);
                img(x, y) = newval;
                err = oldval - newval;

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
} // namespace improcessing
