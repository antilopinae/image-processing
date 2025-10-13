#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <image.hpp>
#include <improcessing.hpp>
#include <png_guard.hpp>
#include <file_closer.hpp>

namespace improcessing {
    auto ReadImage(const std::string &filename) -> std::expected<Image, boost::system::error_code> {
        std::unique_ptr<FILE, details::FileCloser> file(
            std::fopen(filename.c_str(), "rb")
        );

        if (!file) {
            return std::unexpected{
                boost::system::errc::make_error_code(
                    boost::system::errc::no_such_file_or_directory
                )
            };
        }

        PngGuard png(false);

        png_init_io(png.png_ptr, file.get());
        png_read_info(png.png_ptr, png.info_ptr);

        png_uint_32 width = 0;
        png_uint_32 height = 0;
        int bit_depth = 0;
        int color_type = 0;

        png_get_IHDR(png.png_ptr, png.info_ptr, &width, &height,
                     &bit_depth, &color_type, NULL, NULL, NULL);

        if (color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_palette_to_rgb(png.png_ptr);

        if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
            png_set_rgb_to_gray_fixed(png.png_ptr, 1, -1, -1);

        if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_strip_alpha(png.png_ptr);

        if (bit_depth == 16)
            png_set_strip_16(png.png_ptr);

        png_read_update_info(png.png_ptr, png.info_ptr);

        std::vector<uint8_t> data(width * height);
        std::vector<png_bytep> row_ptrs(height);
        for (size_t y = 0; y < height; ++y) {
            row_ptrs[y] = data.data() + y * width;
        }

        png_read_image(png.png_ptr, row_ptrs.data());
        png_read_end(png.png_ptr, png.info_ptr);

        return Image{
            width,
            height,
            std::move(data)
        };
    }

    auto SaveImage(const std::string &filename, const Image &image) -> std::expected<void, boost::system::error_code> {
        if (image.data.size() != image.width * image.height) {
            return std::unexpected{
                boost::system::errc::make_error_code(
                    boost::system::errc::invalid_argument
                )
            };
        }

        std::unique_ptr<FILE, details::FileCloser> file(
            std::fopen(filename.c_str(), "wb")
        );

        if (!file) {
            return std::unexpected{
                boost::system::errc::make_error_code(
                    boost::system::errc::io_error
                )
            };
        }

        PngGuard png(true);

        png_init_io(png.png_ptr, file.get());

        png_set_IHDR(
            png.png_ptr, png.info_ptr,
            image.width, image.height,
            8,
            PNG_COLOR_TYPE_GRAY,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
        );

        png_write_info(png.png_ptr, png.info_ptr);

        std::vector<png_bytep> row_ptrs(image.height);

        for (size_t y = 0; y < image.height; ++y) {
            row_ptrs[y] = const_cast<png_bytep>(
                image.data.data() + y * image.width
            );
        }

        png_write_image(png.png_ptr, row_ptrs.data());
        png_write_end(png.png_ptr, png.info_ptr);

        return {};
    }

    auto MakeCircularGrayscale(Image &image,
                               double radius_fraction) -> std::expected<void, boost::system::error_code> {
        const auto cx = (image.width - 1) / 2.0;
        const auto cy = (image.height - 1) / 2.0;
        const auto r = std::min(image.width, image.height) * radius_fraction;
        const auto edge_width = std::max(1.0, r * 0.02);

        for (int y = 0; y < image.height; ++y) {
            for (int x = 0; x < image.width; ++x) {
                const auto dx = x - cx;
                const auto dy = y - cy;
                const auto d = std::sqrt(dx * dx + dy * dy);

                uint8_t value = 0;

                if (d <= r - edge_width) {
                    const auto t = d / (r - edge_width);
                    auto v = 255.0 * (1.0 - 0.7 * t);

                    v = std::clamp(v, 0.0, 255.0);
                    value = static_cast<uint8_t>(std::round(v));
                } else if (d <= r + edge_width) {
                    const auto t = (d - (r - edge_width)) / (2.0 * edge_width);
                    auto v = 255.0 * (1.0 - t);

                    v = std::clamp(v, 0.0, 255.0);
                    value = static_cast<uint8_t>(std::round(v));
                }

                image.data[y * image.width + x] = value;
            }
        }

        return {};
    }

    auto Blend(const Image &A, const Image &B, const Image &Alpha) -> std::expected<Image, boost::system::error_code> {
        if (A.width > B.width || A.height > B.height ||
            A.width > Alpha.width || A.height > Alpha.height) {
            return std::unexpected{
                boost::system::errc::make_error_code(
                    boost::system::errc::invalid_argument
                )
            };
        }

        const auto totalPixels = static_cast<size_t>(A.width) * static_cast<size_t>(A.height);

        Image out;
        out.width = A.width;
        out.height = A.height;
        out.data.resize(totalPixels);

        for (size_t i = 0; i < totalPixels; ++i) {
            const auto a = A.data[i];
            const auto b = B.data[i];
            const auto alpha = Alpha.data[i];

            const auto blended =
                    static_cast<uint64_t>(alpha) * b +
                    static_cast<uint64_t>(255 - alpha) * a;

            out.data[i] = static_cast<uint8_t>(blended / 255);
        }

        return out;
    }

    auto FloydSteinbergDither(const Image &src, uint8_t n_levels)
        -> std::expected<Image, boost::system::error_code> {
        if (n_levels < 2) {
            return std::unexpected{
                boost::system::errc::make_error_code(boost::system::errc::invalid_argument)
            };
        }

        std::vector<float> buf(src.height * src.width);
        for (size_t i = 0; i < buf.size(); ++i) {
            buf[i] = static_cast<float>(src.data[i]);
        }

        Image out(src.width, src.height);

        const auto levels_minus1 = static_cast<float>(n_levels - 1);
        const auto scale_to_level = levels_minus1 / 255.0f;
        const auto scale_from_level = 255.0f / levels_minus1;

        for (size_t y = 0; y < src.height; ++y) {
            for (size_t x = 0; x < src.width; ++x) {
                const size_t idx = y * src.width + x;
                auto oldv = buf[idx];

                auto level_f = std::round(oldv * scale_to_level);

                if (level_f < 0.0f) level_f = 0.0f;
                if (level_f > levels_minus1) level_f = levels_minus1;

                auto newv = level_f * scale_from_level;
                auto err = oldv - newv;

                uint64_t newbyte = std::lround(newv);

                if (newbyte > 255u) {
                    newbyte = 255u;
                }

                out.data[idx] = static_cast<uint8_t>(newbyte);

                if (x + 1 < src.width) {
                    buf[idx + 1] += err * (7.0f / 16.0f);
                }
                // bottom-left: (x-1, y+1) -> 3/16
                if (x > 0 && (y + 1) < src.height) {
                    buf[idx + src.width - 1] += err * (3.0f / 16.0f);
                }
                // bottom: (x, y+1) -> 5/16
                if ((y + 1) < src.height) {
                    buf[idx + src.width] += err * (5.0f / 16.0f);
                }
                // bottom-right: (x+1, y+1) -> 1/16
                if ((x + 1) < src.width && (src.height + 1) < src.height) {
                    buf[idx + src.width + 1] += err * (1.0f / 16.0f);
                }
            }
        }

        return out;
    }
}
