#include <png.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <image.hpp>

auto SaveGrayPng(const std::string &filename, const Image &image) -> void {
    // if ((int) image.data.size() != image.width * image.height) {
    //     throw std::runtime_error("buffer size mismatch");
    // }
    //
    // FILE *fp = fopen(filename.c_str(), "wb");
    // if (!fp) throw std::runtime_error("Failed to open file for writing: " + filename);
    //
    // png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    // if (!png_ptr) {
    //     fclose(fp);
    //     throw std::runtime_error("Failed to create png_write_struct");
    // }
    //
    // png_infop info_ptr = png_create_info_struct(png_ptr);
    // if (!info_ptr) {
    //     png_destroy_write_struct(&png_ptr, nullptr);
    //     fclose(fp);
    //     throw std::runtime_error("Failed to create png_info_struct");
    // }
    //
    // if (setjmp(png_jmpbuf(png_ptr))) {
    //     png_destroy_write_struct(&png_ptr, &info_ptr);
    //     fclose(fp);
    //     throw std::runtime_error("libpng error during init_io");
    // }
    //
    // png_init_io(png_ptr, fp);
    //
    // // Заголовок: ширина, высота, bit_depth=8, цветовой тип = grayscale
    // png_set_IHDR(png_ptr, info_ptr,
    //              image.width, image.height,
    //              8,
    //              PNG_COLOR_TYPE_GRAY,
    //              PNG_INTERLACE_NONE,
    //              PNG_COMPRESSION_TYPE_DEFAULT,
    //              PNG_FILTER_TYPE_DEFAULT);
    //
    // png_write_info(png_ptr, info_ptr);
    //
    // // Создаём массив указателей на строки
    // std::vector<png_bytep> row_ptrs(image.height);
    // for (int y = 0; y < image.height; ++y) {
    //     row_ptrs[y] = const_cast<png_bytep>(image.data.data() + y * image.width);
    // }
    //
    // png_write_image(png_ptr, row_ptrs.data());
    // png_write_end(png_ptr, nullptr);
    //
    // png_destroy_write_struct(&png_ptr, &info_ptr);
    // fclose(fp);
}

// Генерирует полутоновое круглое изображение (градиент внутри круга)
auto MakeCircularGrayscale(const Image &image, double radius_fraction = 0.45) -> void {
    // auto buf = &image.data;
    // *buf = std::vector<unsigned char>(image.width * image.height, 0);
    //
    // const double cx = (image.width - 1) / 2.0;
    // const double cy = (image.height - 1) / 2.0;
    // const double r = std::min(image.width, image.height) * radius_fraction;
    //
    // // Для плавного края используем небольшую ширину перехода (anti-alias)
    // const double edge_width = std::max(1.0, r * 0.02); // 2% радиуса или минимум 1px
    //
    // for (int y = 0; y < image.height; ++y) {
    //     for (int x = 0; x < image.width; ++x) {
    //         double dx = x - cx;
    //         double dy = y - cy;
    //         double d = std::sqrt(dx * dx + dy * dy);
    //
    //         unsigned char value = 0;
    //         if (d <= r - edge_width) {
    //             // внутри круга (внутреннее тело) — яркость задаём по радиальному градиенту
    //             // пример: яркость уменьшается линейно от центра (255) к краю (64)
    //             double t = d / (r - edge_width); // 0..1
    //             double v = 255.0 * (1.0 - 0.7 * t); // внутри: 255 -> 255*(1-0.7)=76.5
    //             if (v < 0.0) v = 0.0;
    //             if (v > 255.0) v = 255.0;
    //             value = static_cast<unsigned char>(std::round(v));
    //         } else if (d <= r + edge_width) {
    //             // сглаживающий переход (anti-alias): плавное затухание до 0
    //             double t = (d - (r - edge_width)) / (2.0 * edge_width); // 0..1
    //             // взять значение на внутренней границе
    //             double inner_v = 255.0 * (1.0 - 0.7 * ((r - edge_width) / (r - edge_width))); // =255*(1-0.7*1)=76.5
    //             double v = inner_v * (1.0 - t); // линейно гасим
    //             if (v < 0.0) v = 0.0;
    //             if (v > 255.0) v = 255.0;
    //             value = static_cast<unsigned char>(std::round(v));
    //         } else {
    //             // вне круга — 0 (чёрный)
    //             value = 0;
    //         }
    //
    //         buf[y * image.width + x] = value;
    //     }
    // }
}

auto ReadGrayPng(const std::string &filename) -> Image {
    // FILE *fp = fopen(filename.c_str(), "rb");
    // if (!fp) throw std::runtime_error("Failed to open file: " + filename);
    //
    // png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    // if (!png_ptr) {
    //     fclose(fp);
    //     throw std::runtime_error("Failed to create read struct");
    // }
    //
    // png_infop info_ptr = png_create_info_struct(png_ptr);
    // if (!info_ptr) {
    //     png_destroy_read_struct(&png_ptr, nullptr, nullptr);
    //     fclose(fp);
    //     throw std::runtime_error("Failed to create info struct");
    // }
    //
    // if (setjmp(png_jmpbuf(png_ptr))) {
    //     png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
    //     fclose(fp);
    //     throw std::runtime_error("libpng read error");
    // }
    //
    // png_init_io(png_ptr, fp);
    // png_read_info(png_ptr, info_ptr);

    u_int64_t width = 0, height = 0;
    // int bit_depth, color_type;
    // png_get_IHDR(png_ptr, info_ptr, &width, &height,
    //              &bit_depth, &color_type, nullptr, nullptr, nullptr);
    //
    // // Преобразуем всё к 8 бит grayscale
    // if (color_type == PNG_COLOR_TYPE_PALETTE)
    //     png_set_palette_to_rgb(png_ptr);
    // if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
    //     png_set_rgb_to_gray_fixed(png_ptr, 1, -1, -1);
    // if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    //     png_set_strip_alpha(png_ptr);
    // if (bit_depth == 16)
    //     png_set_strip_16(png_ptr);
    //
    // png_read_update_info(png_ptr, info_ptr);

    std::vector<uint8_t> data(width * height);
    // std::vector<png_bytep> row_ptrs(height);
    // for (size_t y = 0; y < height; ++y) {
    //     row_ptrs[y] = data.data() + y * width;
    // }
    //
    // png_read_image(png_ptr, row_ptrs.data());
    // png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
    // fclose(fp);


    return Image{width, height, std::move(data)};
}

auto WriteGrayPng(const std::string &filename, const Image &img) -> void {
    // FILE *fp = fopen(filename.c_str(), "wb");
    // if (!fp) throw std::runtime_error("Failed to open file for writing: " + filename);
    //
    // png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    // png_infop info_ptr = png_create_info_struct(png_ptr);
    // if (!png_ptr || !info_ptr) {
    //     fclose(fp);
    //     throw std::runtime_error("Failed to create png structs");
    // }
    //
    // if (setjmp(png_jmpbuf(png_ptr))) {
    //     png_destroy_write_struct(&png_ptr, &info_ptr);
    //     fclose(fp);
    //     throw std::runtime_error("libpng write error");
    // }
    //
    // png_init_io(png_ptr, fp);
    // png_set_IHDR(png_ptr, info_ptr,
    //              img.width, img.height,
    //              8,
    //              PNG_COLOR_TYPE_GRAY,
    //              PNG_INTERLACE_NONE,
    //              PNG_COMPRESSION_TYPE_DEFAULT,
    //              PNG_FILTER_TYPE_DEFAULT);
    // png_write_info(png_ptr, info_ptr);
    //
    // std::vector<png_bytep> row_ptrs(img.height);
    // for (int y = 0; y < img.height; ++y)
    //     row_ptrs[y] = (png_bytep) (img.data.data() + y * img.width);
    //
    // png_write_image(png_ptr, row_ptrs.data());
    // png_write_end(png_ptr, nullptr);
    // png_destroy_write_struct(&png_ptr, &info_ptr);
    // fclose(fp);
}

auto Blend(const Image &A, const Image &B, const Image &Alpha) -> Image {
    if (A.width != B.width || A.height != B.height ||
        A.width != Alpha.width || A.height != Alpha.height) {
        throw std::runtime_error("Image sizes do not match");
    }

    Image out;
    out.width = A.width;
    out.height = A.height;
    // out.data.resize(out.width * out.height);

    // for (size_t i = 0; i < out.data.size(); ++i) {
    //     unsigned char a = A.data[i];
    //     unsigned char b = B.data[i];
    //     unsigned char alpha = Alpha.data[i];
    //
    //     unsigned int result = (alpha * b + (255 - alpha) * a) / 255;
    //     out.data[i] = static_cast<unsigned char>(result);
    // }

    return out;
}
