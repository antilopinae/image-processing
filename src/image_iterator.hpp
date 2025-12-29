#pragma once

#include <algorithm>
#include <utility>

namespace improcessing {
template<bool IsConst, typename T>
class ImageIterator {
public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;

    explicit ImageIterator(value_type *ptr) : ptr_(ptr) {}

    value_type &operator*() const
    {
        return *ptr_;
    }

    value_type *operator->() const
    {
        return ptr_;
    }

    ImageIterator &operator++()
    {
        for (int i = 0; i < step_; ++i) {
            ++ptr_;
        }
        return *this;
    }

    ImageIterator operator++(int)
    {
        ImageIterator temp = *this;
        for (int i = 0; i < step_; ++i) {
            ++(*this);
        }
        return temp;
    }

    ImageIterator &operator--()
    {
        for (int i = 0; i < step_; ++i) {
            --ptr_;
        }
        return *this;
    }

    ImageIterator operator--(int)
    {
        ImageIterator temp = *this;
        for (int i = 0; i < step_; ++i) {
            --(*this);
        }
        return temp;
    }

    ImageIterator &operator+=(difference_type n)
    {
        ptr_ += n * step_;
        return *this;
    }

    ImageIterator operator+(difference_type n) const
    {
        ImageIterator temp = *this;
        return temp += n * step_;
    }

    ImageIterator &operator-=(difference_type n)
    {
        ptr_ -= n * step_;
        return *this;
    }

    ImageIterator operator-(difference_type n) const
    {
        ImageIterator temp = *this;
        return temp -= n * step_;
    }

    difference_type operator-(const ImageIterator &other) const
    {
        return (ptr_ - other.ptr_) / step_;
    }

    bool operator==(const ImageIterator &other) const
    {
        return ptr_ == other.ptr_;
    }

    bool operator!=(const ImageIterator &other) const
    {
        return !(*this == other);
    }

    bool operator<(const ImageIterator &other) const
    {
        return ptr_ < other.ptr_;
    }

    bool operator>(const ImageIterator &other) const
    {
        return ptr_ > other.ptr_;
    }

    bool operator<=(const ImageIterator &other) const
    {
        return ptr_ <= other.ptr_;
    }

    bool operator>=(const ImageIterator &other) const
    {
        return ptr_ >= other.ptr_;
    }

    template<typename U>
    explicit operator ImageIterator<IsConst, U>() const
    {
        auto v = ImageIterator<IsConst, U>(reinterpret_cast<U *>(ptr_));
        v.set_step(sizeof(U) / sizeof(T));
        static_assert(sizeof(U) > sizeof(T), "Cast Type U must be larger than T");
        return v;
    }

private:
    template<bool B, typename U>
    friend class ImageIterator;

    void set_step(difference_type step) noexcept
    {
        step_ = step;
    }

    value_type *ptr_;
    difference_type step_ = 1;
};
} // namespace improcessing
