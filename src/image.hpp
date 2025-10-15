#pragma once

#include <cinttypes>
#include <cstdint>
#include <memory>
#include <image_vector.hpp>
#include <fmt/base.h>

namespace improcessing {
    template<typename T>
    class ImageBase {
    public:
        using value_type = T;
        using iterator_type = ImVector<value_type>::iterator;
        using const_iterator_type = ImVector<value_type>::const_iterator;
        using vector = ImVector<value_type>;

        explicit ImageBase(size_t width = 0, size_t height = 0)
            : width_(width)
              , height_(height)
              , data_(width * height) {
        }

        ImageBase(size_t width, size_t height, ImVector<value_type> &&data)
            : width_(width)
              , height_(height)
              , data_(std::forward<ImVector<value_type> >(data)) {
        }

        [[nodiscard]] auto get_ptrs() const -> ImVector<value_type *> {
            ImVector<value_type *> row_ptrs(height_);

            for (size_t y = 0; y < height_; ++y) {
                row_ptrs[y] = const_cast<value_type *>(data()) + y * width_;
            }

            return row_ptrs;
        }

        value_type &operator()(size_t x, size_t y) {
            return data_.at(y * width_ + x);
        }

        const value_type &operator()(size_t x, size_t y) const {
            return data_.at(y * width_ + x);
        }

        void resize(size_t new_width, size_t new_height) {
            data_.resize(new_width * new_height);
        }

        iterator_type begin() noexcept {
            return data_.begin();
        }

        iterator_type end() noexcept {
            return data_.end();
        }

        const_iterator_type begin() const noexcept {
            return data_.begin();
        }

        const_iterator_type end() const noexcept {
            return data_.end();
        }

        const_iterator_type cbegin() const noexcept {
            return data_.cbegin();
        }

        const_iterator_type cend() const noexcept {
            return data_.cend();
        }

        [[nodiscard]] size_t size() const noexcept {
            return data_.size();
        }

        [[nodiscard]] bool empty() const noexcept {
            return data_.empty();
        }

        value_type *data() noexcept {
            return data_.data();
        }

        const value_type *data() const noexcept {
            return data_.data();
        }

        [[nodiscard]] size_t width() const noexcept {
            return width_;
        }

        [[nodiscard]] size_t height() const noexcept {
            return height_;
        }

    private:
        ImVector<T> data_;
        size_t width_;
        size_t height_;
    };

    template<typename T = uint8_t>
    class ImageGray : public ImageBase<T> {
    public:
        using ImageBase<T>::ImageBase;
        using ImageBase<T>::value_type;
        using ImageBase<T>::iterator_type;
        using ImageBase<T>::const_iterator_type;
        using ImageBase<T>::vector;
    };

    template<typename T>
    struct PixelRGB {
        union {
            struct {
                T r;
                T g;
                T b;
            };

            T data[3];
        };

        explicit PixelRGB() : r(0), g(0), b(0) {
        }

        PixelRGB(T r_, T g_, T b_) : r(r_), g(g_), b(b_) {
        }
    };

    template<typename T = uint8_t>
    class ImageRGB {
    public:
        enum class Type {
            kRGB,
            kGray
        };

        using ImageGray = ImageBase<T>;

        using gray_type = ImageGray::value_type;
        using vector_gray = ImageGray::vector;
        using iterator_gray_type = ImageGray::iterator_type;
        using const_iterator_gray_type = ImageGray::const_iterator_type;

        using rgb_type = PixelRGB<T>;
        using vector_rgb = ImageBase<rgb_type>;
        using iterator_rgb_type = GenericIterator<false, rgb_type>;
        using const_iterator_rgb_type = GenericIterator<true, rgb_type>;

        explicit ImageRGB(size_t width = 0, size_t height = 0, Type type = Type::kGray) : image_(
            GetWidthFrom(width, type), height) {
        }

        ImageRGB(size_t width, size_t height, vector_gray &&data, Type type = Type::kGray) : image_(
            GetWidthFrom(width, type), height, std::move(data)) {
        }

        [[nodiscard]] auto GetGrayPtrs() const -> ImVector<gray_type *> {
            return image_.get_ptrs();
        }

        [[nodiscard]] auto GetRgbPtrs() const -> ImVector<rgb_type *> {
            auto gray_ptrs = image_.get_ptrs();
            ImVector<rgb_type *> vec(gray_ptrs.size());
            memcpy(vec, gray_ptrs, sizeof(vec));

            return vec;
        }

        auto operator()(size_t x, size_t y) -> gray_type & {
            return image_(x, y);
        }

        auto operator()(size_t x, size_t y) const -> gray_type const & {
            return image_(x, y);
        }

        auto ResizeGray(size_t new_width, size_t new_height) -> void {
            image_.resize(new_width, new_height);
        }

        auto ResizeRgb(size_t new_width, size_t new_height) -> void {
            image_.resize(new_width * kAmountRgb, new_height);
        }

        auto begin() noexcept -> iterator_rgb_type {
            return static_cast<iterator_rgb_type>(image_.begin());
        }

        auto end() noexcept -> iterator_rgb_type {
            return static_cast<iterator_rgb_type>(image_.end());
        }

        auto begin() const noexcept -> const_iterator_rgb_type {
            return static_cast<const_iterator_rgb_type>(image_.begin());
        }

        auto end() const noexcept -> const_iterator_rgb_type {
            return static_cast<const_iterator_rgb_type>(image_.end());
        }

        auto cbegin() const noexcept -> const_iterator_rgb_type {
            return static_cast<const_iterator_rgb_type>(image_.cbegin());
        }

        auto cend() const noexcept -> const_iterator_rgb_type {
            return static_cast<const_iterator_rgb_type>(image_.cend());
        }

        [[nodiscard]] size_t size() const noexcept {
            return image_.size();
        }

        [[nodiscard]] bool Empty() const noexcept {
            return image_.empty();
        }

        gray_type *DataGray() noexcept {
            return image_.data();
        }

        const gray_type *DataGray() const noexcept {
            return image_.data();
        }

        rgb_type *DataRgb() noexcept {
            return static_cast<rgb_type *>(image_.data());
        }

        const rgb_type *DataRgb() const noexcept {
            return static_cast<rgb_type *>(image_.data());
        }

        [[nodiscard]] size_t WidthGray() const noexcept {
            return image_.width();
        }

        [[nodiscard]] size_t HeightGray() const noexcept {
            return image_.height();
        }

        [[nodiscard]] size_t WidthRgb() const noexcept {
            return image_.width() / kAmountRgb;
        }

        [[nodiscard]] size_t HeightRgb() const noexcept {
            return image_.height();
        }

    private:
        static constexpr auto kAmountRgb = sizeof(rgb_type) / sizeof(T);

        constexpr auto GetWidthFrom(size_t from, Type type) const noexcept -> size_t {
            switch (type) {
                case Type::kRGB:
                    return from * kAmountRgb;
                default:
                    return from;
            }
        }

        ImageGray image_;
    };

    using Image = ImageRGB<uint8_t>;
} // namespace improcessing
