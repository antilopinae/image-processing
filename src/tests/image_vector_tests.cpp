#include <gtest/gtest.h>
#include <image_vector.hpp>

using namespace improcessing;

TEST(ImVector, DefaultConstructorCreatesEmptyVector_UnitTest) {
    ImVector<int> v;
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.size(), 0);
    EXPECT_EQ(v.capacity(), 0);
}

TEST(ImVector, ConstructWithSize_UnitTest) {
    ImVector<int> v(5);
    EXPECT_EQ(v.size(), 5);
    EXPECT_EQ(v.capacity(), 5);
}

TEST(ImVector, ConstructWithSizeAndValue_UnitTest) {
    ImVector<int> v(3, 42);
    EXPECT_EQ(v.size(), 3);
    for (auto val: v) {
        EXPECT_EQ(val, 42);
    }
}

TEST(ImVector, InitializerListConstructor_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    EXPECT_EQ(v.size(), 3);
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
    std::cout << v << std::endl;
}

TEST(ImVector, PushBackIncreasesSize_UnitTest) {
    ImVector<int> v;
    v.push_back(5);
    EXPECT_EQ(v.size(), 1);
    EXPECT_EQ(v[0], 5);
}

TEST(ImVector, EmplaceBack_UnitTest) {
    struct Dummy {
        int a;

        explicit Dummy(int x) : a(x) {
        }

        explicit Dummy() = default;
    };
    ImVector<Dummy> v;
    v.emplace_back(123);
    EXPECT_EQ(v.size(), 1);
    EXPECT_EQ(v[0].a, 123);
}

TEST(ImVector, PopBackReducesSize_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    v.pop_back();
    EXPECT_EQ(v.size(), 2);
    EXPECT_EQ(v.back(), 2);
}

TEST(ImVector, ReserveIncreasesCapacity_UnitTest) {
    ImVector<int> v;
    v.reserve(10);
    EXPECT_GE(v.capacity(), 10);
}

TEST(ImVector, ResizeLarger_UnitTest) {
    ImVector<int> v;
    v.resize(5, 7);
    EXPECT_EQ(v.size(), 5);
    for (auto val: v)
        EXPECT_EQ(val, 7);
}

TEST(ImVector, ResizeSmaller_UnitTest) {
    ImVector<int> v = {1, 2, 3, 4, 5};
    v.resize(2);
    EXPECT_EQ(v.size(), 2);
    EXPECT_EQ(v[1], 2);
}

TEST(ImVector, FrontAndBackWork_UnitTest) {
    ImVector<int> v = {10, 20, 30};
    EXPECT_EQ(v.front(), 10);
    EXPECT_EQ(v.back(), 30);
}

TEST(ImVector, OutOfRangeThrows_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    EXPECT_THROW(v.at(3), std::out_of_range);
    EXPECT_THROW(v[3], std::out_of_range);
}

TEST(ImVector, ClearEmptiesVector_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    v.clear();
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.size(), 0);
}

TEST(ImVector, CopyConstructor_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    ImVector<int> copy(v);
    EXPECT_EQ(copy.size(), v.size());
    EXPECT_EQ(copy[0], 1);
}

TEST(ImVector, MoveConstructor_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    ImVector<int> moved(std::move(v));
    EXPECT_EQ(moved.size(), 3);
    EXPECT_TRUE(v.empty());
}

TEST(ImVector, AssignmentOperator_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    ImVector<int> u;
    u = v;
    EXPECT_EQ(u.size(), 3);
    EXPECT_EQ(u[1], 2);
}

TEST(ImVector, MoveAssignmentOperator_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    ImVector<int> u;
    u = std::move(v);
    EXPECT_EQ(u.size(), 3);
    EXPECT_TRUE(v.empty());
}

TEST(ImVectorIterator, ForwardIteration_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    int expected = 1;
    for (auto it = v.begin(); it != v.end(); ++it, ++expected) {
        EXPECT_EQ(*it, expected);
    }
}

TEST(ImVectorIterator, ReverseIteration_UnitTest) {
    ImVector<int> v = {1, 2, 3};
    int expected = 3;
    for (auto it = v.rbegin(); it != v.rend(); ++it, --expected) {
        EXPECT_EQ(*it, expected);
    }
}

TEST(ImVectorIterator, IteratorArithmetic_UnitTest) {
    ImVector<int> v = {10, 20, 30, 40};
    auto it = v.begin();
    it += 2;
    EXPECT_EQ(*it, 30);
    it = it - 1;
    EXPECT_EQ(*it, 20);
    auto diff = v.end() - v.begin();
    EXPECT_EQ(diff, 4);
}

TEST(ImVectorIterator, IteratorComparison_UnitTest) {
    ImVector<int> v = {10, 20, 30};
    auto it1 = v.begin();
    auto it2 = v.begin() + 2;
    EXPECT_TRUE(it1 < it2);
    EXPECT_TRUE(it2 > it1);
    EXPECT_TRUE(it1 != it2);
    EXPECT_FALSE(it1 == it2);
}

TEST(ImVector, Swap_UnitTest) {
    ImVector<int> a = {1, 2};
    ImVector<int> b = {10, 20, 30};
    swap(a, b);
    EXPECT_EQ(a.size(), 3);
    EXPECT_EQ(b.size(), 2);
    EXPECT_EQ(a[0], 10);
    EXPECT_EQ(b[0], 1);
}
