#include <gtest/gtest.h>

TEST(SimpleTest, Addition) {
    auto Add = [](int a, int b) -> int { return a + b; };
    EXPECT_EQ(Add(2, 3), 5);
    EXPECT_EQ(Add(-1, 1), 0);
    EXPECT_NE(Add(2, 2), 5);
}
