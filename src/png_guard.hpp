#pragma once

#include <png.h>

namespace improcessing {
    class PngGuard {
    public:
        png_structp png_ptr = nullptr;
        png_infop info_ptr = nullptr;
        bool is_write = false;

        explicit PngGuard(png_structp p = nullptr, png_infop i = nullptr, bool write = false) noexcept
            : png_ptr(p), info_ptr(i), is_write(write) {
        }

        explicit PngGuard(bool write) : is_write(write) {
            if (is_write) {
                png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
                info_ptr = png_create_info_struct(png_ptr);
            } else {
                png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
                info_ptr = png_create_info_struct(png_ptr);
            }
        }

        ~PngGuard() noexcept {
            if (png_ptr) {
                if (is_write) {
                    png_destroy_write_struct(&png_ptr, &info_ptr);
                } else {
                    png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
                }
            }
        }

        PngGuard(PngGuard &&o) noexcept {
            png_ptr = o.png_ptr;
            info_ptr = o.info_ptr;
            is_write = o.is_write;
            o.png_ptr = nullptr;
            o.info_ptr = nullptr;
        }

        PngGuard &operator=(PngGuard &&o) noexcept {
            if (this != &o) {
                this->~PngGuard();
                png_ptr = o.png_ptr;
                info_ptr = o.info_ptr;
                is_write = o.is_write;
                o.png_ptr = nullptr;
                o.info_ptr = nullptr;
            }
            return *this;
        }

    private:
        PngGuard(const PngGuard &) = delete;

        PngGuard &operator=(const PngGuard &) = delete;
    };
}
