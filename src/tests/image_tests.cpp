#include <gtest/gtest.h>
#include <cstdint>
#include <image.hpp>
#include <numeric>

using namespace improcessing;

TEST(ImageBase, DefaultConstructorCreatesEmptyImage_UnitTest) {
    ImageBase<uint8_t> img;
    EXPECT_EQ(img.width(), 0);
    EXPECT_EQ(img.height(), 0);
    EXPECT_TRUE(img.empty());
    EXPECT_EQ(img.size(), 0);
}

TEST(ImageBase, ConstructorAllocatesCorrectSize_UnitTest) {
    ImageBase<uint8_t> img(4, 5);
    EXPECT_EQ(img.width(), 4);
    EXPECT_EQ(img.height(), 5);
    EXPECT_EQ(img.size(), 20);
    EXPECT_FALSE(img.empty());
}

TEST(ImageBase, ElementAccessAndModification_UnitTest) {
    ImageBase<uint8_t> img(3, 3);
    img(1, 1) = 42;
    EXPECT_EQ(img(1, 1), 42);
}

TEST(ImageBase, ResizeChangesSize_UnitTest) {
    ImageBase<uint8_t> img(2, 2);
    img.resize(5, 5);
    EXPECT_EQ(img.size(), 25);
}

TEST(ImageBase, GetPtrsReturnsRowPointers_UnitTest) {
    ImageBase<uint8_t> img(3, 2);
    img(0, 0) = 11;
    img(0, 1) = 22;

    auto row_ptrs = img.get_ptrs();
    ASSERT_EQ(row_ptrs.size(), 2);
    EXPECT_EQ(row_ptrs[0][0], 11);
    EXPECT_EQ(row_ptrs[1][0], 22);
}

TEST(ImageBase, IteratorsWork_UnitTest) {
    ImageBase<uint8_t> img(2, 2);
    std::iota(img.begin(), img.end(), 1);
    int sum = 0;
    for (auto v: img) sum += v;
    EXPECT_EQ(sum, 10); // 1+2+3+4
}

TEST(ImageRGB, DefaultConstructor_UnitTest) {
    Image img;
    EXPECT_TRUE(img.Empty());
    EXPECT_EQ(img.WidthGray(), 0);
    EXPECT_EQ(img.HeightGray(), 0);
}

TEST(ImageRGB, GrayImageInitialization_UnitTest) {
    Image img(4, 3, Image::Type::kGray);
    EXPECT_EQ(img.WidthGray(), 4);
    EXPECT_EQ(img.HeightGray(), 3);
    EXPECT_EQ(img.size(), 12);
}

TEST(ImageRGB, RGBImageInitialization_UnitTest) {
    Image img(2, 2, Image::Type::kRGB);
    EXPECT_EQ(img.WidthRgb(), 2);
    EXPECT_EQ(img.HeightRgb(), 2);
    EXPECT_EQ(img.WidthGray(), 2 * 3); // 3 channels
    EXPECT_EQ(img.size(), 2 * 2 * 3);
}

TEST(ImageRGB, ElementWriteAndReadGray_UnitTest) {
    Image img(2, 2, Image::Type::kGray);
    img(1, 1) = 77;
    EXPECT_EQ(img(1, 1), 77);
}

TEST(ImageRGB, ResizeGray_UnitTest) {
    Image img(2, 2);
    img.ResizeGray(4, 3);
    EXPECT_EQ(img.WidthGray(), 4);
    EXPECT_EQ(img.HeightGray(), 3);
    EXPECT_EQ(img.size(), 12);
}

TEST(ImageRGB, ResizeRgb_UnitTest) {
    Image img(2, 2, Image::Type::kRGB);
    img.ResizeRgb(5, 4);
    EXPECT_EQ(img.WidthRgb(), 5);
    EXPECT_EQ(img.HeightRgb(), 4);
    EXPECT_EQ(img.size(), 5 * 4 * 3);
}

TEST(ImageRGB, GetGrayPtrsReturnsValidData_UnitTest) {
    Image img(2, 2, Image::Type::kGray);
    img(0, 0) = 10;
    auto gray_ptrs = img.GetGrayPtrs();
    EXPECT_EQ(gray_ptrs[0][0], 10);
}

TEST(ImageRGB, GetRGBPtrsReturnsValidData_UnitTest) {
    Image img(3, 1, Image::Type::kGray);
    img(0, 0) = 1;
    img(1, 0) = 2;
    img(2, 0) = 3;
    auto rgb_ptrs = img.GetRgbPtrs();
    auto pixel = PixelRGB<uint8_t>{1, 2, 3};
    EXPECT_EQ(rgb_ptrs[0][0], pixel);
}

TEST(ImageRGB, DataRgbPointerCastsCorrectly_UnitTest) {
    Image img(3, 2, Image::Type::kRGB);
    img(0, 0) = 1;
    img(1, 0) = 2;
    img(2, 0) = 3;

    auto *rgb_data = img.DataRgb();

    EXPECT_EQ(rgb_data[0].r, 1);
    EXPECT_EQ(rgb_data[0].g, 2);
    EXPECT_EQ(rgb_data[0].b, 3);
}

TEST(ImageRGB, EmptyReturnsTrueWhenNoData_UnitTest) {
    Image img;
    EXPECT_TRUE(img.Empty());
}

TEST(ImageRGB, Iterators_UnitTest) {
    Image img(3, 2, Image::Type::kRGB);
    img(0, 0) = 1;
    img(1, 0) = 2;
    img(2, 0) = 3;

    img(0, 1) = 1;
    img(1, 1) = 2;
    img(2, 1) = 3;

    int i = 0;
    for (auto v: img) {
        EXPECT_EQ(v.r, 1);
        EXPECT_EQ(v.g, 2);
        EXPECT_EQ(v.b, 3);
        ++i;
    }

    EXPECT_EQ(i, 2);

    i = 0;
    for (const auto v: img) {
        EXPECT_EQ(v.r, 1);
        EXPECT_EQ(v.g, 2);
        EXPECT_EQ(v.b, 3);
        ++i;
    }

    EXPECT_EQ(i, 2);
}
