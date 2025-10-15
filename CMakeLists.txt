cmake_minimum_required(VERSION 3.26.0 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CCache OPTIONAL)
include(NoInSourceBuilds OPTIONAL)

option(USE_VCPKG_TOOLCHAIN "Use vcpkg toolchain file" ON)

if (USE_VCPKG_TOOLCHAIN)
    # set(VCPKG_OVERLAY_PORTS "${CMAKE_SOURCE_DIR}/vcpkg_ports")

    if (NOT DEFINED CMAKE_TOOLCHAIN_FILE)
        set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE FILEPATH "Vcpkg toolchain file")
    endif ()
endif ()

if (WIN32 AND MINGW)
    set(VCPKG_TARGET_TRIPLET "x64-mingw-static" CACHE STRING "vcpkg triplet" FORCE)
elseif (WIN32 AND MSVC)
    set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "vcpkg triplet" FORCE)
endif()

# project setup
project(image-processing VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

# dependencies
# find_package(PkgConfig REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(fmt REQUIRED)
find_package(argparse CONFIG REQUIRED)
find_package(PNG REQUIRED)

# pkg_check_modules(LIBBPF REQUIRED IMPORTED_TARGET libbpf)

add_subdirectory(src)
add_subdirectory(src/tests)