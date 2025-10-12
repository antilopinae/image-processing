#pragma once

#include <vector>
#include <functional>
#include <memory>

struct Image {
    Image(u_int64_t width = 0, u_int64_t height = 0) : width(width), height(height) {
    }

    Image(u_int64_t width, u_int64_t height, std::vector<uint8_t> &&data) : width(width), height(height),
                                                                            data(std::forward<std::vector<uint8_t> >(
                                                                                data)) {
    }

    u_int64_t width;
    u_int64_t height;
    std::vector<uint8_t> data;
};
