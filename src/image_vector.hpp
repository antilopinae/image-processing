#pragma once

#include <memory>
#include <algorithm>
#include <stdexcept>
#include <concepts>
#include <utility>
#include <functional>
#include <image_iterator.hpp>

namespace improcessing {
    template<typename T>
    class ImVector {
    public:
        // remove ref and others from T - fuck
        using value_type = T;

        using iterator = ImageIterator<false, value_type>;
        using const_iterator = ImageIterator<true, value_type>;
        using reverse_iterator = std::reverse_iterator<ImageIterator<false, value_type> >;
        using const_reverse_iterator = std::reverse_iterator<ImageIterator<true, value_type> >;

        explicit ImVector() : size_(0), capacity_(0), data_(nullptr) {
        }

        explicit ImVector(size_t count) : size_(count), capacity_(count) {
            if (count > 0) {
                data_ = std::make_unique<value_type[]>(capacity_);
            }
        }

        ImVector(size_t count, size_t capacity, value_type *data) noexcept : size_(count), capacity_(capacity),
                                                                             data_(data) {
        }

        ImVector(size_t count, const value_type &value) : size_(count), capacity_(count) {
            if (count > 0) {
                data_ = std::make_unique<value_type[]>(capacity_);
                for (size_t i = 0; i < size_; ++i) {
                    *(data_.get() + i) = T(value);
                }
            }
        }

        ImVector(std::initializer_list<value_type> init) : size_(init.size()), capacity_(init.size()) {
            if (init.size() > 0) {
                data_ = std::make_unique<value_type[]>(capacity_);
                std::copy(init.begin(), init.end(), data_.get());
            }
        }

        ImVector(const ImVector &other) : size_(other.size_), capacity_(other.capacity_) {
            if (other.size_ > 0) {
                data_ = std::make_unique<value_type[]>(capacity_);
                std::copy(other.data_.get(), other.data_.get() + other.size_, data_.get());
            }
        }

        ImVector(ImVector &&other) noexcept
            : size_(std::exchange(other.size_, 0)),
              capacity_(std::exchange(other.capacity_, 0)),
              data_(std::exchange(other.data_, nullptr)) {
        }

        auto operator=(const ImVector &other) -> ImVector & {
            if (this != &other) {
                if (capacity_ < other.size_) {
                    data_ = std::make_unique<value_type[]>(other.capacity_);
                    capacity_ = other.capacity_;
                }
                size_ = other.size_;
                std::copy(other.data_.get(), other.data_.get() + other.size_, data_.get());
            }
            return *this;
        }

        auto operator=(ImVector &&other) noexcept -> ImVector & {
            if (this != &other) {
                data_ = std::exchange(other.data_, nullptr);
                size_ = std::exchange(other.size_, 0);
                capacity_ = std::exchange(other.capacity_, 0);
            }
            return *this;
        }

        auto operator=(std::initializer_list<value_type> ilist) -> ImVector & {
            if (capacity_ < ilist.size()) {
                data_ = std::make_unique<value_type[]>(ilist.size());
                capacity_ = ilist.size();
            }
            size_ = ilist.size();
            std::copy(ilist.begin(), ilist.end(), data_.get());
            return *this;
        }

        auto operator[](size_t index) -> value_type & {
            if (index >= size_) throw std::out_of_range("ImVector::[]: index out of range");
            return data_[index];
        }

        auto operator[](size_t index) const -> const value_type & {
            if (index >= size_) throw std::out_of_range("ImVector::[]: index out of range");
            return data_[index];
        }

        auto at(size_t index) -> value_type & {
            if (index >= size_) throw std::out_of_range("ImVector::at: index out of range");
            return data_[index];
        }

        auto at(size_t index) const -> const value_type & {
            if (index >= size_) throw std::out_of_range("ImVector::at: index out of range");
            return data_[index];
        }

        auto front() -> value_type & {
            return data_[0];
        }

        auto front() const -> const value_type & {
            return data_[0];
        }

        auto back() -> value_type & {
            return data_[size_ - 1];
        }

        auto back() const -> const value_type & {
            return data_[size_ - 1];
        }

        auto data() noexcept -> value_type * {
            return data_.get();
        }

        auto data() const noexcept -> const value_type * {
            return data_.get();
        }

        auto begin() noexcept -> iterator {
            return iterator(data_.get());
        }

        auto begin() const noexcept -> const_iterator {
            return const_iterator(data_.get());
        }

        auto cbegin() const noexcept -> const_iterator {
            return const_iterator(data_.get());
        }

        auto end() noexcept -> iterator {
            return iterator(data_.get() + size_);
        }

        auto end() const noexcept -> const_iterator {
            return const_iterator(data_.get() + size_);
        }

        auto cend() const noexcept -> const_iterator {
            return const_iterator(data_.get() + size_);
        }

        auto rend() noexcept -> reverse_iterator {
            return reverse_iterator(begin());
        }

        auto rbegin() noexcept -> reverse_iterator {
            return reverse_iterator(end());
        }

        auto rbegin() const noexcept -> const_reverse_iterator {
            return const_reverse_iterator(cend());
        }

        auto crbegin() const noexcept -> const_reverse_iterator {
            return const_reverse_iterator(cend());
        }

        auto rend() const noexcept -> const_reverse_iterator {
            return const_reverse_iterator(cbegin());
        }

        auto crend() const noexcept -> const_reverse_iterator {
            return const_reverse_iterator(cbegin());
        }

        auto empty() const noexcept -> bool {
            return size_ == 0;
        }

        auto size() const noexcept -> size_t {
            return size_;
        }

        auto capacity() const noexcept -> size_t {
            return capacity_;
        }

        auto reserve(size_t new_capacity) -> void {
            if (new_capacity > capacity_) {
                auto new_data = std::make_unique<value_type[]>(new_capacity);

                if (data_) {
                    for (size_t i = 0; i < size_; ++i) {
                        *(new_data.get() + i) = value_type(std::move(data_[i]));
                    }
                }

                data_ = std::move(new_data);
                capacity_ = new_capacity;
            }
        }

        auto clear() noexcept -> void {
            for (size_t i = 0; i < size_; ++i) {
                data_[i].~T();
            }
            size_ = 0;
        }

        auto push_back(const value_type &value) -> void {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = value_type(value);
            size_++;
        }

        auto push_back(value_type &&value) -> void {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = value_type(std::move(value));
            size_++;
        }

        template<typename... Args>
        auto emplace_back(Args &&... args) -> value_type & {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = value_type(std::forward<Args>(args)...);
            size_++;
            return data_[size_ - 1];
        }

        auto pop_back() -> void {
            data_[size_ - 1].~T();
            size_--;
        }

        auto resize(size_t count) -> void {
            resize(count, value_type{});
        }

        auto resize(size_t count, const value_type &value) -> void {
            if (count < size_) {
                for (size_t i = count; i < size_; ++i) {
                    data_[i].~T();
                }
                size_ = count;
            } else if (count > size_) {
                if (count > capacity_) {
                    reserve(count);
                }
                for (size_t i = size_; i < count; ++i) {
                    *(data_.get() + i) = value_type(value);
                }
                size_ = count;
            }
        }

        auto swap(ImVector &other) noexcept (std::is_nothrow_swappable_v<std::unique_ptr<value_type[]> >) -> void {
            using std::swap;
            swap(data_, other.data_);
            swap(size_, other.size_);
            swap(capacity_, other.capacity_);
        }

        auto release() noexcept -> value_type * {
            size_ = 0;
            capacity_ = 0;
            return data_.release();
        }

        auto reset() noexcept -> void {
            size_ = 0;
            capacity_ = 0;
            data_.reset();
        }

    private:
        size_t size_;
        size_t capacity_;

        std::unique_ptr<value_type []> data_;
    };

    template<typename T>
    void swap(ImVector<T> &lhs, ImVector<T> &rhs) noexcept(noexcept(lhs.swap(rhs))) {
        lhs.swap(rhs);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &os, const ImVector<T> &vec) {
        // os << "ImVector [size=" << vec.size() << ", capacity=" << vec.capacity() << "] { ";
        // for (const auto &val: vec) {
        //     os << val << " ";
        // }
        // os << "}";
        return os;
    }
}
