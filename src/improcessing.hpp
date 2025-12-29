#pragma once

#include <expected>
#include <string>

#include <image.hpp>
#include <point.hpp>

#include <boost/system/errc.hpp>

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
auto SaveImage(const std::string &filename, const Image &image, Image::Type type = Image::Type::kGray)
    -> std::expected<void, boost::system::error_code>;

/*!
 * @brief Prints circular grayscale at the Image
 * @param image Image in which circle will be printed
 * @param radius_fraction
 * @return nothing or error code in failure
 */
auto MakeCircularGrayscale(Image &image, double radius_fraction = 0.45)
    -> std::expected<void, boost::system::error_code>;

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
auto DrawLine(Image &img, Point2D start, Point2D end, Pixel color) -> std::expected<void, boost::system::error_code>;

/**
 * @brief Fills the circle with the specified color
 * @param img Drawing image
 * @param center Coordinates of the center (Point)
 * @param radius Radius
 * @param color Pixel color
 */
auto FillCircle(Image &img, Point center, double radius, Pixel color) -> void;

enum class LineCap {
    kButt,
    kSquare,
    kRound
};

/*!
 * @brief Draw thick line in Image
 * @param img Image for drawing
 * @param start Coordinates of the beginning of the line
 * @param end Coordinates of the end of the line
 * @param thickness Thickness
 * @param color Color of the line
 * @param cap Line cap
 * @return nothing or error code in failure
 */
auto DrawThickLine(Image &img, Point2D start, Point2D end, double thickness, Pixel color, LineCap cap = LineCap::kButt)
    -> void;

/*!
 * @brief Draw polygon in Image
 * @param img Image for drawing
 * @param poly Coordinates of the vertexes of polygon
 * @param color Color of the edges
 * @return nothing or error code in failure
 */
auto DrawPolygonEdges(Image &img, const std::vector<Point2D> &poly, Pixel color)
    -> std::expected<void, boost::system::error_code>;

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
 *   Edges are considered inside
 * @param img Image for drawing
 * @param poly Coordinates of the vertexes of polygon
 * @param color Color of the filling
 * @return nothing or error code in failure
 */
auto FillPolygonEvenOdd(Image &img, const std::vector<Point2D> &poly, Pixel color) -> void;

/*!
 * @brief Filling the polygon using the non-zero-winding rule for determining whether a pixel belongs to a polygon
 *   Edges are considered inside
 * @param img Image for drawing
 * @param poly Coordinates of the vertexes of polygon
 * @param color Color of the filling
 * @return nothing or error code in failure
 */
auto FillPolygonNonZero(Image &img, const std::vector<Point2D> &poly, Pixel color) -> void;

/*!
 * @brief The Cyrus-Beck algorithm
 * @param p0 coordinates of the beginning of segment line
 * @param p1 coordinates of the end of segment line
 * @param poly Coordinates of the vertexes of polygon that clips line
 * @note param p0 and p1 will be changed to borders of clipped segment
 * @return true if at least part of the segment remains in the polygon, false — if completely outside
 */
auto CyrusBeckClipSegment(Point &p0, Point &p1, const std::vector<Point> &poly) -> bool;

/*!
 * @brief Sutherland-Hodgman polygon clipping
 * Clips an arbitrary subject polygon by a convex clip polygon
 */
auto ClipPolygonSutherlandHodgman(const std::vector<Point> &subjectPoly, const std::vector<Point> &clipPoly)
    -> std::vector<Point>;

/*!
 * @brief The Bezier algorithm for drawing curve
 *   Step count is calculated based on control point distance
 * @param p0, p1, p2, p3 Points to make bezier cubic curve
 * @return vector of steps+1 points of the constructed curve
 */
auto BezierCubicCurve(Point p0, Point p1, Point p2, Point p3) -> std::vector<Point>;

enum class ProjectionType {
    kParallel,
    kPerspective
};

/*!
 * @brief Renders a 3D parallelepiped with rotation, projection and hidden edges
 * @param img Target image
 * @param center Center of the object in 3D space
 * @param size Dimensions
 * @param rotationAxis Axis of rotation
 * @param angle Rotation angle in radians
 * @param type Projection type
 * @param k Distance to the center of projection
 */
auto RenderParallelepiped(Image &img,
                          Point3 center,
                          Point3 size,
                          Point3 rotationAxis,
                          double angle,
                          ProjectionType type,
                          double k) -> void;

/*!
 * @brief Extracts the external contour of a self-intersecting polygon
 *      Uses rasterization (NZW) and Moore-Neighbor tracing
 * @param poly Original polygon
 * @param width Canvas width for processing
 * @param height Canvas height for processing
 * @return A vector of points representing the external contour
 */
auto ExtractExteriorContour(const std::vector<Point2D> &poly, size_t width, size_t height) -> std::vector<Point2D>;

/*!
 * @brief Approximation of a circle arc using Cubic Bezier curves
 * @param center Center of the circle
 * @param radius Radius
 * @param angle_start_deg Start angle in degrees
 * @param angle_end_deg End angle in degrees
 * @return Vector of points forming the arc
 */
auto MakeCircleArc(Point center, double radius, double angle_start_deg, double angle_end_deg) -> std::vector<Point>;

/*!
 * @brief Color quantization using K-Means clustering
 * @param img Image to process
 * @param k Number of colors
 */
auto ColorQuantizationKMeans(Image &img, int k) -> std::expected<void, boost::system::error_code>;

enum class PerspectiveType {
    kTwoPoint,
    kThreePoint
};

struct SceneObject {
    Point3 center;
    Point3 size;
    Point3 rotation_axis;
    double rotation_angle;
    std::vector<Pixel> face_colors;
    PerspectiveType type;
};

auto RenderLab5Scene(Image &img, const SceneObject &obj1, const SceneObject &obj2, double k) -> void;
} // namespace improcessing
