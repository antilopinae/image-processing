find_package(GTest REQUIRED)

enable_testing()

add_executable(${PROJECT_NAME}-tests
    image_vector_tests.cpp
    image_tests.cpp
    floyd_tests.cpp
)

target_include_directories(${PROJECT_NAME}-tests PUBLIC
    ${CMAKE_SOURCE_DIR}/src/tests
)

target_link_libraries(${PROJECT_NAME}-tests PUBLIC
    ${PROJECT_NAME}-core
    GTest::gtest
    GTest::gtest_main
)
