#pragma once

#include <cinttypes>
#include <cstdint>
#include <memory>
#include <vector>

namespace improcessing {
    class Image {
    public:
        using container_type = std::vector<uint8_t>;
        using size_type = container_type::size_type;
        using iterator = container_type::iterator;
        using const_iterator = container_type::const_iterator;

        explicit Image(size_type width = 0, size_type height = 0)
            : width_(width)
              , height_(height)
              , data_(width * height, 0) {
        }

        Image(size_type width, size_type height, std::vector<uint8_t> &&data)
            : width_(width)
              , height_(height)
              , data_(std::forward<std::vector<uint8_t> >(data)) {
        }

        uint8_t &operator()(size_type x, size_type y) {
            return data_.at(y * width_ + x);
        }

        const uint8_t &operator()(size_type x, size_type y) const {
            return data_.at(y * width_ + x);
        }

        void resize(size_type new_width, size_type new_height) {
            data_.resize(new_width * new_height, 0);
        }

        [[nodiscard]] auto get_ptrs() const -> std::vector<uint8_t *> {
            std::vector<uint8_t *> row_ptrs(height_);

            for (size_type y = 0; y < height_; ++y) {
                row_ptrs[y] = const_cast<uint8_t *>(data_.data()) + y * width_;
            }

            return std::move(row_ptrs);
        }

        iterator begin() noexcept { return data_.begin(); }
        iterator end() noexcept { return data_.end(); }
        const_iterator begin() const noexcept { return data_.begin(); }
        const_iterator end() const noexcept { return data_.end(); }
        const_iterator cbegin() const noexcept { return data_.cbegin(); }
        const_iterator cend() const noexcept { return data_.cend(); }

        [[nodiscard]] size_type size() const noexcept { return data_.size(); }
        [[nodiscard]] bool empty() const noexcept { return data_.empty(); }

        uint8_t *data() noexcept { return data_.data(); }
        const uint8_t *data() const noexcept { return data_.data(); }

        [[nodiscard]] size_type width() const noexcept { return width_; }
        [[nodiscard]] size_type height() const noexcept { return height_; }

    private:
        size_type width_;
        size_type height_;
        std::vector<uint8_t> data_;
    };
} // namespace improcessing
