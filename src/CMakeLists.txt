add_library(${PROJECT_NAME}-core STATIC
    improcessing.cpp
)

target_include_directories(${PROJECT_NAME}-core PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(${PROJECT_NAME}-core PUBLIC
    Boost::system
    fmt::fmt-header-only
    argparse::argparse
    PNG::PNG
)

# main executable
add_executable(${PROJECT_NAME}
    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_link_options(${PROJECT_NAME} PRIVATE -static)
if(UNIX)
    target_link_options(${PROJECT_NAME} PRIVATE -static-libgcc -static-libstdc++)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${PROJECT_NAME}-core
)
