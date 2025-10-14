#pragma once

#include <cinttypes>
#include <cstdint>
#include <memory>
#include <vector>

namespace improcessing {
    class Image {
    public:
        explicit Image(std::uint64_t width = 0, std::uint64_t height = 0)
            : width(width)
              , height(height)
              , data(width * height, static_cast<uint8_t>(0)) {
        }

        Image(std::uint64_t width, std::uint64_t height, std::vector<uint8_t> &&data)
            : width(width)
              , height(height)
              , data(std::forward<std::vector<uint8_t> >(data)) {
        }

        std::uint64_t width;
        std::uint64_t height;
        std::vector<uint8_t> data;
    };
} // namespace improcessing
