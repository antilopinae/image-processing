#pragma once

#include <boost/system/errc.hpp>
#include <expected>
#include <string>
#include <image.hpp>

namespace improcessing {
    /*!
     * @brief Reads Image from file
     * @param filename Name of file
     * @return Image at success
     */
    auto ReadImage(const std::string &filename) -> std::expected<Image, boost::system::error_code>;

    /*!
     * @brief Save Image to file
     * @param filename Name of file
     * @return nothing or error code in failure
     */
    auto SaveImage(const std::string &filename, const Image &image) -> std::expected<void, boost::system::error_code>;

    /*!
     * @brief Prints circular grayscale at the Image
     * @param image Image in which circle will be printed
     * @param radius_fraction
     * @return nothing or error code in failure
     */
    auto MakeCircularGrayscale(Image &image,
                               double radius_fraction = 0.45) -> std::expected<void, boost::system::error_code>;

    /*!
     * @brief Blending two images, uses third as alpha channel
     * @param A Image for blending
     * @param B Image for blending
     * @param Alpha Image which will be used as alpha channel
     * @return new Image or error code in failure
     */
    auto Blend(const Image &A, const Image &B, const Image &Alpha) -> std::expected<Image, boost::system::error_code>;

    /*!
     * @brief The Floyd-Stenberg error scattering algorithm
     * @param src Image for blending
     * @param n_levels amount of bpp to convert
     * @return new scatter Image or error code in failure
     */
    auto FloydSteinbergDither(const Image &src, uint8_t n_levels) -> std::expected<Image, boost::system::error_code>;
} // namespace improcessing
