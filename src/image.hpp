#pragma once

#include <vector>
#include <memory>

namespace improcessing {
    class Image {
    public:
        explicit Image(u_int64_t width = 0, u_int64_t height = 0) : width(width), height(height),
                                                                    data(width * height, static_cast<uint8_t>(0)) {
        }

        Image(u_int64_t width, u_int64_t height, std::vector<uint8_t> &&data) : width(width), height(height),
            data(std::forward<std::vector<uint8_t> >(
                data)) {
        }

        u_int64_t width;
        u_int64_t height;
        std::vector<uint8_t> data;
    };
}
