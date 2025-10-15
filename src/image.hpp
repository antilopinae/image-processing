#pragma once

#include <cinttypes>
#include <cstdint>
#include <memory>
#include <vector>

namespace improcessing {
    template<typename T, typename S = std::vector<T>::size_type>
    class ImageBase {
    public:
        using value_type = T;
        using size_type = S;
        using iterator_type = std::vector<value_type>::iterator;
        using const_iterator_type = std::vector<value_type>::const_iterator;

        explicit ImageBase(size_type width = 0, size_type height = 0)
            : width_(width)
              , height_(height)
              , data_(width * height) {
        }


        ImageBase(size_type width, size_type height, std::vector<value_type> &&data)
            : width_(width)
              , height_(height)
              , data_(std::forward<std::vector<value_type> >(data)) {
        }

        [[nodiscard]] auto get_ptrs() const -> std::vector<value_type *> {
            std::vector<value_type *> row_ptrs(height_);

            for (size_type y = 0; y < height_; ++y) {
                row_ptrs[y] = const_cast<value_type *>(data()) + y * width_;
            }

            return std::move(row_ptrs);
        }

        value_type &operator()(size_type x, size_type y) {
            return data_.at(y * width_ + x);
        }

        const value_type &operator()(size_type x, size_type y) const {
            return data_.at(y * width_ + x);
        }

        void resize(size_type new_width, size_type new_height) {
            data_.resize(new_width * new_height);
        }

        iterator_type begin() noexcept { return data_.begin(); }
        iterator_type end() noexcept { return data_.end(); }

        const_iterator_type begin() const noexcept { return data_.begin(); }
        const_iterator_type end() const noexcept { return data_.end(); }

        const_iterator_type cbegin() const noexcept { return data_.cbegin(); }
        const_iterator_type cend() const noexcept { return data_.cend(); }

        [[nodiscard]] size_type size() const noexcept { return data_.size(); }
        [[nodiscard]] bool empty() const noexcept { return data_.empty(); }

        value_type *data() noexcept { return data_.data(); }
        const value_type *data() const noexcept { return data_.data(); }

        [[nodiscard]] size_type width() const noexcept { return width_; }
        [[nodiscard]] size_type height() const noexcept { return height_; }

    private:
        explicit ImageBase(value_type *ptr, size_type width = 0, size_type height = 0)
            : width_(width)
              , height_(height) {
            data_.resize(width * height);
            auto c = data_.data();
            const auto &a = data_.data();
            auto &b = const_cast<value_type *&>(a);
            b = ptr;
            delete[] c;
        }

        friend class Image;

        std::vector<value_type> data_;
        size_type width_;
        size_type height_;
    };

    class Image255 : public ImageBase<uint8_t> {
    public:
        using ImageBase::size_type;
        using ImageBase::ImageBase;
    };

    template<typename T = uint8_t>
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

    class Image255RGB : public ImageBase<PixelRGB<uint8_t> > {
    public:
        using ImageBase::ImageBase;
    };

    class Image {
    public:
        enum class Type {
            kRGB,
            kGray
        };

        using gray_type = Image255::value_type;
        using size_type = Image255::size_type;
        using rgb_type = Image255RGB::value_type;
        using iterator_type = Image255RGB::iterator_type;
        using const_iterator_type = Image255RGB::const_iterator_type;

        static constexpr auto kAmountRgb = sizeof(rgb_type) / sizeof(gray_type);

        explicit Image(size_type width = 0, size_type height = 0, Type type = Type::kGray) {
            if (type == Type::kRGB) {
                width = width * kAmountRgb;
                height = height * kAmountRgb;
            }

            new(&image_.gray_) Image255(width, height);
        }

        Image(size_type width, size_type height, std::vector<gray_type> &&data) {
            new(&image_.gray_) Image255(width, height, std::move(data));
        }

        Image(size_type width, size_type height, std::vector<rgb_type> &&data) {
            new(&image_.gray_) Image255((gray_type *) data.data(), width * kAmountRgb,
                                        height * kAmountRgb);
        }

        Image(Image &&other) noexcept {
            new(&image_.gray_) Image255(std::move(other.image_.gray_));
        }

        Image &operator=(Image &&other) noexcept {
            if (this != &other) {
                image_.gray_.~Image255();
                new(&image_.gray_) Image255(std::move(other.image_.gray_));
            }

            return *this;
        }

        ~Image() {
            image_.gray_.~Image255();
        }

        [[nodiscard]] auto get_ptrs() const -> std::vector<gray_type *> {
            return image_.gray_.get_ptrs();
        }

        rgb_type &operator()(size_type x, size_type y) {
            return image_.rgb_(x, y);
        }

        const rgb_type &operator()(size_type x, size_type y) const {
            return image_.rgb_(x, y);
        }

        void resize(size_type new_width, size_type new_height, Type type = Type::kGray) {
            if (type == Type::kRGB) {
                image_.rgb_.resize(new_width, new_height);
            } else {
                image_.gray_.resize(new_width, new_height);
            }
        }

        iterator_type begin() noexcept { return image_.rgb_.begin(); }
        iterator_type end() noexcept { return image_.rgb_.end(); }

        const_iterator_type begin() const noexcept { return image_.rgb_.begin(); }
        const_iterator_type end() const noexcept { return image_.rgb_.end(); }

        const_iterator_type cbegin() const noexcept { return image_.rgb_.cbegin(); }
        const_iterator_type cend() const noexcept { return image_.rgb_.cend(); }

        [[nodiscard]] size_type size() const noexcept { return image_.gray_.size(); }

        [[nodiscard]] bool empty(Type type = Type::kGray) const noexcept {
            if (type == Type::kRGB) {
                return image_.rgb_.empty();
            }

            return image_.gray_.empty();
        }

        gray_type *data() noexcept { return image_.gray_.data(); }
        const gray_type *data() const noexcept { return image_.gray_.data(); }

        [[nodiscard]] size_type width(Type type = Type::kGray) const noexcept {
            if (type == Type::kRGB) {
                return image_.rgb_.width_;
            }

            return image_.gray_.width_;
        }

        [[nodiscard]] size_type height(Type type = Type::kGray) const noexcept {
            if (type == Type::kRGB) {
                return image_.rgb_.height_;
            }

            return image_.gray_.height_;
        }

    private:
        union Storage {
            Image255RGB rgb_;
            Image255 gray_;

            Storage() {
            }

            ~Storage() {
            }
        } image_;
    };
} // namespace improcessing
