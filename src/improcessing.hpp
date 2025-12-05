#pragma once

#include <boost/system/errc.hpp>
#include <expected>
#include <string>
#include <image.hpp>
#include <point.hpp>

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
    * @param type Type of the Image (if RGB 3 bits per pixel were used)
    * @return nothing or error code in failure
    */
    auto SaveImage(const std::string &filename, const Image &image,
                   Image::Type type = Image::Type::kGray) -> std::expected<void, boost::system::error_code>;

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
    * @param src Image for dithering
    * @param n_levels amount of bpp to convert
    * @return nothing or error code in failure
    */
    auto FloydSteinbergDither(Image &src, uint8_t n_levels) -> std::expected<void, boost::system::error_code>;

    /*!
    * @brief Draw line in Image
    * @param img Image for drawing
    * @param start Coordinates of the beginning of the line
    * @param end Coordinates of the end of the line
    * @param color Color of the line
    * @return nothing or error code in failure
    */
    auto DrawLine(Image &img, Point2D start, Point2D end,
                  Pixel color) -> std::expected<void, boost::system::error_code>;

    /*!
    * @brief Draw polygon in Image
    * @param img Image for drawing
    * @param poly Coordinates of the vertexes of polygon
    * @param color Color of the edges
    * @return nothing or error code in failure
    */
    auto DrawPolygonEdges(Image &img, const std::vector<Point2D> &poly,
                          Pixel color) -> std::expected<void, boost::system::error_code>;

    /*!
    * @brief Determine if the polygon simple or complex (i.e. with self-intersections)
    * @param poly Coordinates of the vertexes of polygon
    * @return true if polygon simple
    */
    auto IsSimplePolygon(const std::vector<Point2D> &poly) -> bool;

    /*!
    * @brief Determine if the polygon convex
    * @param poly Coordinates of the vertexes of polygon
    * @return true if polygon is convex or error code in failure
    */
    auto IsConvexPolygon(const std::vector<Point2D> &poly) -> bool;

    /*!
    * @brief Filling the polygon using the even-odd rule for determining whether a pixel belongs to a polygon
    * @param img Image for drawing
    * @param poly Coordinates of the vertexes of polygon
    * @param color Color of the filling
    * @return nothing or error code in failure
    */
    auto FillPolygonEvenOdd(Image &img, const std::vector<Point2D> &poly,
                            Pixel color) -> void;

    /*!
    * @brief Filling the polygon using the non-zero-winding rule for determining whether a pixel belongs to a polygon
    * @param img Image for drawing
    * @param poly Coordinates of the vertexes of polygon
    * @param color Color of the filling
    * @return nothing or error code in failure
    */
    auto FillPolygonNonZero(Image &img, const std::vector<Point2D> &poly,
                            Pixel color) -> void;
} // namespace improcessing
