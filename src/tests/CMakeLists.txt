find_package(GTest REQUIRED)

enable_testing()

add_executable(${PROJECT_NAME}-tests
        simple_tests.cpp
        floyd_tests.cpp
        #        image_tests.cpp # redundant
)

target_link_libraries(${PROJECT_NAME}-tests PRIVATE GTest::gtest GTest::gtest_main)
