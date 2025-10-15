#pragma once

#include <memory>
#include <algorithm>
#include <stdexcept>
#include <concepts>
#include <utility>
#include <functional>

namespace improcessing {
    template<bool IsConst, typename T>
    class GenericIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;

        explicit GenericIterator(value_type *ptr) : ptr_(ptr) {
        }

        value_type &operator*() const {
            return *ptr_;
        }

        value_type *operator->() const {
            return ptr_;
        }

        GenericIterator &operator++() {
            ++ptr_;
            return *this;
        }

        GenericIterator operator++(int) {
            GenericIterator temp = *this;
            ++(*this);
            return temp;
        }

        GenericIterator &operator--() {
            --ptr_;
            return *this;
        }

        GenericIterator operator--(int) {
            GenericIterator temp = *this;
            --(*this);
            return temp;
        }

        GenericIterator &operator+=(difference_type n) {
            ptr_ += n;
            return *this;
        }

        GenericIterator operator+(difference_type n) const {
            GenericIterator temp = *this;
            return temp += n;
        }

        GenericIterator &operator-=(difference_type n) {
            ptr_ -= n;
            return *this;
        }

        GenericIterator operator-(difference_type n) const {
            GenericIterator temp = *this;
            return temp -= n;
        }

        difference_type operator-(const GenericIterator &other) const { return ptr_ - other.ptr_; }

        bool operator==(const GenericIterator &other) const { return ptr_ == other.ptr_; }
        bool operator!=(const GenericIterator &other) const { return !(*this == other); }
        bool operator<(const GenericIterator &other) const { return ptr_ < other.ptr_; }
        bool operator>(const GenericIterator &other) const { return ptr_ > other.ptr_; }
        bool operator<=(const GenericIterator &other) const { return ptr_ <= other.ptr_; }
        bool operator>=(const GenericIterator &other) const { return ptr_ >= other.ptr_; }

        template<typename G>
        explicit operator G() const {
            return reinterpret_cast<G *>(const_cast<G *>(ptr_));
        }

    private:
        value_type *ptr_;
    };

    template<typename T>
    class ImVector {
    public:
        using value_type = T;

        using iterator = GenericIterator<false, value_type>;
        using const_iterator = GenericIterator<true, value_type>;
        using reverse_iterator = std::reverse_iterator<GenericIterator<false, value_type> >;
        using const_reverse_iterator = std::reverse_iterator<GenericIterator<true, value_type> >;

        explicit ImVector() : size_(0), capacity_(0), data_(nullptr) {
        }

        explicit ImVector(size_t count) : size_(count), capacity_(count) {
            if (count > 0) {
                data_ = std::make_unique<T[]>(capacity_);
            }
        }

        ImVector(size_t count, const value_type &value) : size_(count), capacity_(count) {
            if (count > 0) {
                data_ = std::make_unique<T[]>(capacity_);
                for (size_t i = 0; i < size_; ++i) {
                    *(data_.get() + i) = T(value);
                }
            }
        }

        ImVector(std::initializer_list<T> init) : size_(init.size()), capacity_(init.size()) {
            if (init.size() > 0) {
                data_ = std::make_unique<T[]>(capacity_);
                std::copy(init.begin(), init.end(), data_.get());
            }
        }

        ImVector(const ImVector &other) : size_(other.size_), capacity_(other.capacity_) {
            if (other.size_ > 0) {
                data_ = std::make_unique<T[]>(capacity_);
                std::copy(other.data_.get(), other.data_.get() + other.size_, data_.get());
            }
        }

        ImVector(ImVector &&other) noexcept
            : size_(std::exchange(other.size_, 0)),
              capacity_(std::exchange(other.capacity_, 0)),
              data_(std::exchange(other.data_, nullptr)) {
        }

        ImVector &operator=(const ImVector &other) {
            if (this != &other) {
                if (capacity_ < other.size_) {
                    data_ = std::make_unique<T[]>(other.capacity_);
                    capacity_ = other.capacity_;
                }
                size_ = other.size_;
                std::copy(other.data_.get(), other.data_.get() + other.size_, data_.get());
            }
            return *this;
        }

        ImVector &operator=(ImVector &&other) noexcept {
            if (this != &other) {
                data_ = std::exchange(other.data_, nullptr);
                size_ = std::exchange(other.size_, 0);
                capacity_ = std::exchange(other.capacity_, 0);
            }
            return *this;
        }

        ImVector &operator=(std::initializer_list<T> ilist) {
            if (capacity_ < ilist.size()) {
                data_ = std::make_unique<T[]>(ilist.size());
                capacity_ = ilist.size();
            }
            size_ = ilist.size();
            std::copy(ilist.begin(), ilist.end(), data_.get());
            return *this;
        }

        value_type &operator[](size_t index) {
            if (index >= size_) throw std::out_of_range("ImVector::[]: index out of range");
            return data_[index];
        }

        const value_type &operator[](size_t index) const {
            if (index >= size_) throw std::out_of_range("ImVector::[]: index out of range");
            return data_[index];
        }

        value_type &at(size_t index) {
            if (index >= size_) throw std::out_of_range("ImVector::at: index out of range");
            return data_[index];
        }

        const value_type &at(size_t index) const {
            if (index >= size_) throw std::out_of_range("ImVector::at: index out of range");
            return data_[index];
        }

        value_type &front() {
            return data_[0];
        }

        const value_type &front() const {
            return data_[0];
        }

        value_type &back() {
            return data_[size_ - 1];
        }

        const value_type &back() const {
            return data_[size_ - 1];
        }

        value_type *data() noexcept {
            return data_.get();
        }

        const value_type *data() const noexcept {
            return data_.get();
        }

        iterator begin() noexcept {
            return iterator(data_.get());
        }

        const_iterator begin() const noexcept {
            return const_iterator(data_.get());
        }

        const_iterator cbegin() const noexcept {
            return const_iterator(data_.get());
        }

        iterator end() noexcept {
            return iterator(data_.get() + size_);
        }

        const_iterator end() const noexcept {
            return const_iterator(data_.get() + size_);
        }

        const_iterator cend() const noexcept {
            return const_iterator(data_.get() + size_);
        }

        reverse_iterator rbegin() noexcept {
            return reverse_iterator(end());
        }

        const_reverse_iterator rbegin() const noexcept {
            return const_reverse_iterator(cend());
        }

        const_reverse_iterator crbegin() const noexcept {
            return const_reverse_iterator(cend());
        }

        reverse_iterator rend() noexcept {
            return reverse_iterator(begin());
        }

        const_reverse_iterator rend() const noexcept {
            return const_reverse_iterator(cbegin());
        }

        const_reverse_iterator crend() const noexcept {
            return const_reverse_iterator(cbegin());
        }

        bool empty() const noexcept {
            return size_ == 0;
        }

        size_t size() const noexcept {
            return size_;
        }

        size_t capacity() const noexcept {
            return capacity_;
        }

        void reserve(size_t new_capacity) {
            if (new_capacity > capacity_) {
                std::unique_ptr<T[]> new_data = std::make_unique<T[]>(new_capacity);

                if (data_) {
                    for (size_t i = 0; i < size_; ++i) {
                        *(new_data.get() + i) = T(std::move(data_[i]));
                    }
                }

                data_ = std::move(new_data);
                capacity_ = new_capacity;
            }
        }

        void clear() noexcept {
            for (size_t i = 0; i < size_; ++i) {
                data_[i].~T();
            }
            size_ = 0;
        }

        void push_back(const value_type &value) {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = T(value);
            size_++;
        }

        void push_back(T &&value) {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = T(std::move(value));
            size_++;
        }

        template<typename... Args>
        value_type &emplace_back(Args &&... args) {
            if (size_ == capacity_) {
                reserve(capacity_ == 0 ? 1 : capacity_ * 2);
            }
            *(data_.get() + size_) = T(std::forward<Args>(args)...);
            size_++;
            return data_[size_ - 1];
        }

        void pop_back() {
            data_[size_ - 1].~T();
            size_--;
        }

        void resize(size_t count) {
            resize(count, T{});
        }

        void resize(size_t count, const value_type &value) {
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
                    *(data_.get() + i) = T(value);
                }
                size_ = count;
            }
        }

        void swap(ImVector &other) noexcept (std::is_nothrow_swappable_v<std::unique_ptr<T[]> >) {
            using std::swap;
            swap(data_, other.data_);
            swap(size_, other.size_);
            swap(capacity_, other.capacity_);
        }

    private:
        size_t size_;
        size_t capacity_;

        std::unique_ptr<T []> data_;
    };

    template<typename T>
    void swap(ImVector<T> &lhs, ImVector<T> &rhs) noexcept(noexcept(lhs.swap(rhs))) {
        lhs.swap(rhs);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &os, const ImVector<T> &vec) {
        os << "ImVector [size=" << vec.size() << ", capacity=" << vec.capacity() << "] { ";
        for (const auto &val: vec) {
            os << val << " ";
        }
        os << "}";
        return os;
    }
}
