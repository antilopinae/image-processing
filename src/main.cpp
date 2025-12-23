#include <cassert>
#include <commands.hpp>
#include <expected>
#include <filesystem>
#include <image.hpp>
#include <improcessing.hpp>
#include <iostream>

#include <argparse/argparse.hpp>
#include <boost/system/error_code.hpp>
#include <fmt/core.h>

auto AddInputArgument(argparse::ArgumentParser &parser, std::string &input,
                      const std::string &name = "--input") -> void {
    parser.add_argument(name).required().store_into(input).help("specify the input file");
}

auto AddOutputArgument(argparse::ArgumentParser &parser, std::string &output) -> void {
    parser.add_argument("-o", "--output").required().store_into(output).help("specify the output file");
}

auto DoSmthWithImage(const std::string &im_name, const std::function<void(improcessing::Image &)> &payload,
                     int width = 600,
                     int height = 300) -> std::expected<void, boost::system::error_code> {
    improcessing::Image img(width, height, improcessing::Image::Type::kRGB);
    payload(img);

    auto save = SaveImage(im_name, img, improcessing::Image::Type::kRGB);
    if (!save) {
        return std::unexpected{save.error()};
    }

    fmt::print("Saved an image: {}\n", (std::filesystem::current_path() / im_name).string());

    return {};
}

int main(int argc, char *argv[]) {
    std::string input_first{};
    std::string input_second{};
    std::string input_third{};
    std::string output{};

    argparse::ArgumentParser parser("image-processing", "0.1.0", argparse::default_arguments::help);
    parser.add_description("Image processing program");

    argparse::ArgumentParser command_lab1("lab1");
    command_lab1.add_description("[circle-gray | blend]");

    argparse::ArgumentParser command_lab1_circle_gray("circle-gray");
    command_lab1_circle_gray.add_description("prints circle gray");
    AddOutputArgument(command_lab1_circle_gray, output);

    argparse::ArgumentParser command_lab1_blending("blend");
    command_lab1_blending.add_description("blending two images");
    AddOutputArgument(command_lab1_blending, output);
    AddInputArgument(command_lab1_blending, input_first, "--image-first");
    AddInputArgument(command_lab1_blending, input_second, "--image-second");
    AddInputArgument(command_lab1_blending, input_third, "--image-alpha");

    command_lab1.add_subparser(command_lab1_circle_gray);
    command_lab1.add_subparser(command_lab1_blending);
    parser.add_subparser(command_lab1);

    uint8_t n_levels;
    argparse::ArgumentParser command_lab2("lab2");
    command_lab2.add_description("starts Floyd-Stenberg algorithm");
    command_lab2.add_argument("-n", "--n-levels").required().store_into(n_levels).help("set level for scattering");
    AddInputArgument(command_lab2, input_first);
    AddOutputArgument(command_lab2, output);

    parser.add_subparser(command_lab2);

    argparse::ArgumentParser command_lab3("lab3");
    command_lab3.add_description("starts Laboratory 3");

    parser.add_subparser(command_lab3);

    argparse::ArgumentParser command_lab4("lab4");
    command_lab4.add_description("starts Laboratory 4");

    parser.add_subparser(command_lab4);

    argparse::ArgumentParser command_lab5("lab5");
    command_lab5.add_description("3D Projections and Animation");

    auto k_param = 200.0;
    auto frames = 36;

    command_lab5.add_argument("-k", "--perspective-k").scan<'g', double>().store_into(k_param).default_value(200.0).
            help("Perspective center Z coordinate");

    command_lab5.add_argument("--axis")
            .nargs(3)
            .scan<'g', double>()
            .default_value(std::vector<double>{1.0, 1.0, 0.0})
            .help("Rotation axis x y z");

    command_lab5.add_argument("-n", "--frames").scan<'i', int>().store_into(frames).default_value(36).help(
        "Number of frames for animation");

    command_lab5.add_argument("-o", "--output-prefix").required().store_into(output).help(
        "Output filename prefix (e.g. 'frame')");

    parser.add_subparser(command_lab5);

    argparse::ArgumentParser command_hw1("hw1");
    command_hw1.add_description("Extract external contour of self-intersecting polygon");

    AddOutputArgument(command_hw1, output);

    parser.add_subparser(command_hw1);

    argparse::ArgumentParser command_hw2("hw2");
    command_hw2.add_description("Draw circular arc using Bezier curves");

    double arc_cx = 300, arc_cy = 300, arc_r = 200;
    double arc_start = 0, arc_end = 270;

    AddOutputArgument(command_hw2, output);

    command_hw2.add_argument("--center").nargs(2).scan<'g', double>().default_value(std::vector<double>{300, 300}).
            help("Center X Y");

    command_hw2.add_argument("-r", "--radius").scan<'g', double>().store_into(arc_r).default_value(200.0);

    command_hw2.add_argument("-s", "--start").scan<'g', double>().store_into(arc_start).default_value(0.0).help(
        "Start angle (deg)");

    command_hw2.add_argument("-e", "--end").scan<'g', double>().store_into(arc_end).default_value(270.0).help(
        "End angle (deg)");

    parser.add_subparser(command_hw2);

    argparse::ArgumentParser command_hw3("hw3");
    command_hw3.add_description("K-Means Color Quantization");

    int k_colors = 8;

    AddInputArgument(command_hw3, input_first);
    AddOutputArgument(command_hw3, output);

    command_hw3.add_argument("-k").scan<'i', int>().store_into(k_colors).default_value(8).help("Number of colors");

    parser.add_subparser(command_hw3);

    try {
        parser.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    auto initialize_cmd = [&] -> std::optional<Command> {
        if (parser.is_subcommand_used(command_lab1)) {
            if (command_lab1.is_subcommand_used(command_lab1_circle_gray)) {
                return Laboratory1GrayCircle{
                    .width = 512, .height = 512, .radius_fraction = 0.45, .output_filename = output
                };
            } else if (command_lab1.is_subcommand_used(command_lab1_blending)) {
                return Laboratory1Blend{
                    .first_filename = input_first, .second_filename = input_second, .alpha_filename = input_third,
                    .output_filename = output
                };
            }
        } else if (parser.is_subcommand_used(command_lab2)) {
            return Laboratory2{.input_filename = input_first, .output_filename = output, .n_levels = n_levels};
        } else if (parser.is_subcommand_used(command_lab3)) {
            return Laboratory3{};
        } else if (parser.is_subcommand_used(command_lab4)) {
            return Laboratory4{};
        } else if (parser.is_subcommand_used(command_lab5)) {
            auto axis = command_lab5.get<std::vector<double> >("--axis");

            return Laboratory5{
                .output_prefix = output,
                .k = k_param,
                .axis_x = axis[0], .axis_y = axis[1], .axis_z = axis[2],
                .frames = frames
            };
        } else if (parser.is_subcommand_used(command_hw1)) {
            return Homework1{.output_filename = output};
        } else if (parser.is_subcommand_used(command_hw2)) {
            auto c = command_hw2.get<std::vector<double> >("--center");

            return Homework2{
                .output_filename = output,
                .center_x = c[0], .center_y = c[1], .radius = arc_r,
                .angle_start = arc_start, .angle_end = arc_end
            };
        } else if (parser.is_subcommand_used(command_hw3)) {
            return Homework3{.input_filename = input_first, .output_filename = output, .k_colors = k_colors};
        }

        return std::nullopt;
    };

    std::optional<Command> cmd = initialize_cmd();

    using namespace improcessing;

    const auto v = overloaded{
        [](const Laboratory1GrayCircle &c) -> std::expected<void, boost::system::error_code> {
            auto image = Image{c.width, c.height};

            auto grayscale = MakeCircularGrayscale(image, c.radius_fraction);
            if (!grayscale) {
                return std::unexpected{grayscale.error()};
            }

            auto save = SaveImage(c.output_filename, image);
            if (!save) {
                return std::unexpected{save.error()};
            }

            fmt::print(
                "Saved grayscale png image: {}, with width {} and height {}, with radius fraction: {}\n",
                (std::filesystem::current_path() / c.output_filename).string(),
                c.width,
                c.height,
                c.radius_fraction);

            return {};
        },
        [](const Laboratory1Blend &c) -> std::expected<void, boost::system::error_code> {
            auto A = ReadImage(c.first_filename);
            if (!A) {
                return std::unexpected{A.error()};
            }

            auto B = ReadImage(c.second_filename);
            if (!B) {
                return std::unexpected{B.error()};
            }

            auto Alpha = ReadImage(c.alpha_filename);
            if (!Alpha) {
                return std::unexpected{Alpha.error()};
            }

            auto blended = Blend(A.value(), B.value(), Alpha.value());
            if (!blended) {
                return std::unexpected{blended.error()};
            }

            auto result = SaveImage(c.output_filename, blended.value());
            if (!result) {
                return std::unexpected{result.error()};
            }

            fmt::print("Saved blended png image: {}\n", (std::filesystem::current_path() / c.output_filename).string());

            return {};
        },
        [](const Laboratory2 &c) -> std::expected<void, boost::system::error_code> {
            auto image = ReadImage(c.input_filename);
            if (!image) {
                return std::unexpected{image.error()};
            }

            auto result = FloydSteinbergDither(image.value(), c.n_levels);
            if (!result) {
                return std::unexpected{result.error()};
            }

            auto save = SaveImage(c.output_filename, image.value());
            if (!save) {
                return std::unexpected{save.error()};
            }

            fmt::print(
                "Saved scatter image: {}, with n_levels: {}\n",
                (std::filesystem::current_path() / c.output_filename).string(), c.n_levels);

            return {};
        },
        [](const Laboratory3 &c) -> std::expected<void, boost::system::error_code> {
            std::vector<Point2D> poly = {
                {0, 0}, {0, 20}, {20, 20}, {20, 0}
            };

            int width = 600, height = 300;

            for (auto &p: poly) {
                p.y += height / 2;
                p.x += width / 2;
            }

            constexpr auto image_with_edges_name = "PolygonEdges.png";
            constexpr auto image_with_filled_polys_name = "PolygonFilled.png";

            Image img(width, height, Image::Type::kRGB);

            std::vector<Point2D> right = poly;
            std::vector<Point2D> left = poly;

            for (auto &p: left) { p.x -= width / 4; }
            for (auto &p: right) { p.x += width / 4; }

            if (auto draw_left_polygon_edges = DrawPolygonEdges(img, left, {255, 230, 100}); !draw_left_polygon_edges) {
                return std::unexpected{draw_left_polygon_edges.error()};
            }

            if (auto draw_right_polygon_edges = DrawPolygonEdges(img, right, {255, 230, 100}); !
                draw_right_polygon_edges) {
                return std::unexpected{draw_right_polygon_edges.error()};
            }

            auto save = SaveImage(image_with_edges_name, img, Image::Type::kRGB);
            if (!save) {
                return std::unexpected{save.error()};
            }

            fmt::print(
                "Saved an image with drawn polygon edges: {}\n",
                (std::filesystem::current_path() / image_with_edges_name).string());

            FillPolygonEvenOdd(img, left, {200, 220, 255});
            FillPolygonNonZero(img, right, {200, 220, 255});

            if (auto draw_left_polygon_edges = DrawPolygonEdges(img, left, {255, 230, 100}); !draw_left_polygon_edges) {
                return std::unexpected{draw_left_polygon_edges.error()};
            }

            if (auto draw_right_polygon_edges = DrawPolygonEdges(img, right, {255, 230, 100}); !
                draw_right_polygon_edges) {
                return std::unexpected{draw_right_polygon_edges.error()};
            }

            save = SaveImage(image_with_filled_polys_name, img, Image::Type::kRGB);
            if (!save) {
                return std::unexpected{save.error()};
            }

            fmt::print(
                "Saved an image with filled polygons: {}\n",
                (std::filesystem::current_path() / image_with_filled_polys_name).string());

            auto simple = IsSimplePolygon(poly);
            auto convex = IsConvexPolygon(poly);

            fmt::print(
                "Determine polygon as {} and {}\n",
                simple ? "simple" : "complex (self-intersecting)", convex ? "convex" : "non-convex");

            {
                int w = 500, h = 500;
                Image imgEO(w, h, Image::Type::kRGB);
                Image imgNZW(w, h, Image::Type::kRGB);

                std::vector<Point2D> star = {
                    {250, 50}, {310, 450}, {50, 200}, {450, 200}, {190, 450}
                };

                FillPolygonEvenOdd(imgEO, star, {100, 150, 255});
                FillPolygonNonZero(imgNZW, star, {255, 100, 100});

                SaveImage("lab3_star_EO.png", imgEO, Image::Type::kRGB);
                SaveImage("lab3_star_NZW.png", imgNZW, Image::Type::kRGB);
                fmt::print("Saved EO vs NZW comparison (lab3_star_*.png)\n");
            }

            {
                std::vector<Point2D> convex = {{10, 10}, {50, 10}, {50, 50}, {10, 50}};
                std::vector<Point2D> concave = {{10, 10}, {50, 10}, {30, 30}, {50, 50}, {10, 50}};
                std::vector<Point2D> complex = {{10, 10}, {50, 50}, {50, 10}, {10, 50}};

                assert(IsConvexPolygon(convex));
                assert(!IsConvexPolygon(concave));
                assert(!IsConvexPolygon(complex));

                fmt::print("Square is convex: {}\n", IsConvexPolygon(convex)); // true
                fmt::print("Arrow is convex: {}\n", IsConvexPolygon(concave)); // false
                fmt::print("Bowtie is simple: {}\n", IsSimplePolygon(complex)); // false
            }

            auto res = DoSmthWithImage("lab3_thick_lines.png", [](Image &img) {
                DrawThickLine(img, {50, 50}, {200, 50}, 10, {255, 0, 0}, LineCap::kButt);
                DrawThickLine(img, {50, 100}, {200, 100}, 10, {0, 255, 0}, LineCap::kSquare);
                DrawThickLine(img, {50, 150}, {200, 150}, 10, {0, 0, 255}, LineCap::kRound);
            }, 300, 250);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab3_thick_lines_demo.png", [](Image &img) {
                img.ResizeRgb(600, 600);

                DrawThickLine(img, {50, 50}, {550, 50}, 20, {255, 100, 100}, LineCap::kButt);
                DrawThickLine(img, {50, 100}, {550, 100}, 20, {100, 255, 100}, LineCap::kSquare);
                DrawThickLine(img, {50, 150}, {550, 150}, 20, {100, 100, 255}, LineCap::kRound);

                DrawThickLine(img, {100, 200}, {100, 550}, 30, {200, 200, 0}, LineCap::kButt);
                DrawThickLine(img, {200, 200}, {200, 550}, 30, {0, 200, 200}, LineCap::kRound);

                DrawThickLine(img, {300, 250}, {550, 500}, 50, {255, 255, 255}, LineCap::kRound);
                DrawThickLine(img, {300, 500}, {550, 250}, 5, {255, 0, 255}, LineCap::kSquare);

                DrawLine(img, {50, 50}, {550, 50}, {0, 0, 0});
                DrawLine(img, {50, 100}, {550, 100}, {0, 0, 0});
                DrawLine(img, {50, 150}, {550, 150}, {0, 0, 0});
            }, 600, 600);

            return {};
        },
        [](const Laboratory4 &c) -> std::expected<void, boost::system::error_code> {
            auto res = DoSmthWithImage("segment_cross_square.png", [](Image &img) -> void {
                std::vector<Point> square = {{5, 0}, {5, 10}, {15, 10}, {15, 0}};
                Point A(0, 5), B(20, 5);

                Point a = A, b = B;
                bool inside = CyrusBeckClipSegment(a, b, square);

                assert(inside && "The segment not cross square!");

                std::ignore = DrawLine(img, A.ToPoint2D(), B.ToPoint2D(), {200, 220, 255});
                std::ignore = DrawLine(img, a.ToPoint2D(), b.ToPoint2D(), {100, 30, 150});

                std::vector<Point2D> square2d;
                square2d.reserve(square.size());

                std::ranges::transform(square.begin(), square.end(), std::back_inserter(square2d), [](auto &p) {
                    return p.ToPoint2D();
                });

                std::ignore = DrawPolygonEdges(img, square2d, {255, 230, 100});
            }, 30, 15);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("segment_inside_square.png", [](Image &img) -> void {
                std::vector<Point> square = {{5, 0}, {5, 10}, {15, 10}, {15, 0}};
                Point A(5, 5), B(10, 10);

                Point a = A, b = B;
                bool inside = CyrusBeckClipSegment(a, b, square);

                assert(inside && "The segment not cross square!");

                std::ignore = DrawLine(img, A.ToPoint2D(), B.ToPoint2D(), {200, 220, 255});
                std::ignore = DrawLine(img, a.ToPoint2D(), b.ToPoint2D(), {100, 30, 150});

                std::vector<Point2D> square2d;
                square2d.reserve(square.size());

                std::ranges::transform(square.begin(), square.end(), std::back_inserter(square2d), [](auto &p) {
                    return p.ToPoint2D();
                });

                std::ignore = DrawPolygonEdges(img, square2d, {255, 230, 100});
            }, 30, 15);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("segment_outside_square.png", [](Image &img) -> void {
                std::vector<Point> square = {{5, 5}, {5, 15}, {15, 15}, {15, 5}};
                Point A(0, 0), B(4, 4);

                Point a = A, b = B;
                bool inside = CyrusBeckClipSegment(a, b, square);

                assert(!inside && "The segment cross square!");

                std::ignore = DrawLine(img, A.ToPoint2D(), B.ToPoint2D(), {200, 220, 255});

                std::vector<Point2D> square2d;
                square2d.reserve(square.size());

                std::ranges::transform(square.begin(), square.end(), std::back_inserter(square2d), [](auto &p) {
                    return p.ToPoint2D();
                });

                std::ignore = DrawPolygonEdges(img, square2d, {255, 230, 100});
            }, 30, 20);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("bezier_symmetric.png", [](Image &img) -> void {
                Point p0(10, 10);
                Point p1(20, 40);
                Point p2(40, 40);
                Point p3(50, 10);

                auto curve = BezierCubicCurve(p0, p1, p2, p3);

                assert(curve.front().Equal(p0));
                assert(curve.back().Equal(p3));

                for (int i = 0; i <= 30; i++) {
                    auto &c1 = curve[i];
                    auto &c2 = curve[60 - i];

                    assert(fabs((c1.x + c2.x)/2 - 30.0) < 0.5);
                }

                for (auto &c: curve) {
                    const auto &p = c.ToPoint2D();
                    img.GetRGBPixel(p.x, p.y) = {200, 255, 200};
                }
            }, 60, 50);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("bezier_linearity.png", [](Image &img) -> void {
                Point p0(0, 0), p1(5, 5), p2(10, 10), p3(15, 15);

                auto curve = BezierCubicCurve(p0, p1, p2, p3);

                for (size_t i = 0; i < curve.size(); i++) {
                    double t = double(i) / 20.0;

                    Point expected = p0 * (1 - t) + p3 * t;
                    assert(curve[i].Equal(expected, 1e-6));
                }

                for (auto &c: curve) {
                    const auto &p = c.ToPoint2D();
                    img.GetRGBPixel(p.x, p.y) = {200, 255, 200};
                }
            }, 20, 20);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("bezier_basic.png", [](Image &img) -> void {
                Point p0(2, 2);
                Point p1(5, 15);
                Point p2(15, 15);
                Point p3(18, 2);

                auto curve = BezierCubicCurve(p0, p1, p2, p3);

                assert(curve.size() == steps + 1);

                assert(curve.front().Equal(p0));
                assert(curve.back().Equal(p3));

                std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {255, 200, 100});
                std::ignore = DrawLine(img, p1.ToPoint2D(), p2.ToPoint2D(), {255, 200, 100});
                std::ignore = DrawLine(img, p2.ToPoint2D(), p3.ToPoint2D(), {255, 200, 100});

                for (auto &c: curve) {
                    const auto &p = c.ToPoint2D();
                    img.GetRGBPixel(p.x, p.y) = {200, 255, 200};
                }
            }, 20, 20);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab4_bezier_adaptive.png", [](Image &img)-> void {
                img.ResizeRgb(600, 600);
                Point p0(50, 550), p1(50, 50), p2(550, 50), p3(550, 550);

                auto curve = BezierCubicCurve(p0, p1, p2, p3);

                for (const auto &p: curve) {
                    auto p2d = p.ToPoint2D();
                    if (p2d.x < 600 && p2d.y < 600)
                        img.GetRGBPixel(p2d.x, p2d.y) = {0, 255, 0};
                }
            });

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab4_cyrus_beck.png", [](Image &img)-> void {
                img.ResizeRgb(500, 500);
                // ccw
                std::vector<Point> clipPoly = {{100, 300}, {300, 300}, {300, 100}, {100, 100}};

                std::vector<Point2D> clip2d;
                for (auto &p: clipPoly) clip2d.push_back(p.ToPoint2D());
                std::ignore = DrawPolygonEdges(img, clip2d, {255, 255, 0});

                std::vector<std::pair<Point, Point> > lines = {
                    {{50, 200}, {350, 200}},
                    {{150, 150}, {250, 250}},
                    {{50, 50}, {400, 400}},
                    {{50, 50}, {80, 80}}
                };

                for (auto [p0, p1]: lines) {
                    std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {50, 50, 50});
                    if (CyrusBeckClipSegment(p0, p1, clipPoly)) {
                        std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {255, 0, 0});
                    }
                }
            });

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab4_sutherland_hodgman.png", [](Image &img)-> void {
                img.ResizeRgb(600, 600);
                std::vector<Point> subject = {{100, 100}, {250, 50}, {400, 100}, {400, 400}, {100, 400}};
                std::vector<Point> clip = {{200, 0}, {500, 250}, {200, 500}, {0, 250}};

                auto result = ClipPolygonSutherlandHodgman(subject, clip);

                std::vector<Point2D> sub2d, clip2d, res2d;
                for (auto &p: subject) sub2d.push_back(p.ToPoint2D());
                for (auto &p: clip) clip2d.push_back(p.ToPoint2D());
                for (auto &p: result) res2d.push_back(p.ToPoint2D());

                std::ignore = DrawPolygonEdges(img, sub2d, {100, 100, 100});
                std::ignore = DrawPolygonEdges(img, clip2d, {255, 255, 0});

                if (!res2d.empty()) {
                    FillPolygonEvenOdd(img, res2d, {0, 255, 0});
                    std::ignore = DrawPolygonEdges(img, res2d, {255, 255, 255});
                }
                fmt::print("Sutherland-Hodgman clipped polygon size: {}\n", result.size());
            });

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab4_degenerate_test.png", [](Image &img) {
                std::vector<Point> degeneratePoly = {{10, 10}, {50, 10}, {30, 10}};
                Point a(0, 0), b(100, 100);
                bool clipped = CyrusBeckClipSegment(a, b, degeneratePoly);
                assert(!clipped && "Should not clip against a line-polygon");
            });

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("lab4_complex_convex_clipping.png", [](Image &img) {
                img.ResizeRgb(600, 600);
                std::vector<Point> clipPoly = {
                    {100, 300}, {150, 100}, {450, 100},
                    {500, 300}, {450, 500}, {150, 500}
                };

                std::vector<Point2D> clip2d;
                for (auto &p: clipPoly) clip2d.push_back(p.ToPoint2D());
                std::ignore = DrawPolygonEdges(img, clip2d, {255, 255, 0});

                std::vector<std::pair<Point, Point> > lines = {
                    {{50, 300}, {550, 300}},
                    {{300, 50}, {300, 550}},
                    {{100, 100}, {500, 500}},
                    {{150, 100}, {450, 500}},
                    {{50, 50}, {150, 100}},
                    {{450, 100}, {550, 100}}
                };

                for (auto [p0, p1]: lines) {
                    std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {60, 60, 60});

                    Point a = p0, b = p1;
                    if (CyrusBeckClipSegment(a, b, clipPoly)) {
                        std::ignore = DrawLine(img, a.ToPoint2D(), b.ToPoint2D(), {0, 255, 0});
                    }
                }
            }, 600, 600);

            res = DoSmthWithImage("lab4_star_sutherland_hodgman.png", [](Image &img) {
                img.ResizeRgb(600, 600);
                std::vector<Point> star = {
                    {300, 50}, {350, 200}, {500, 200}, {380, 300},
                    {450, 450}, {300, 350}, {150, 450}, {220, 300},
                    {100, 200}, {250, 200}
                };

                std::vector<Point> clip = {
                    {300, 100}, {500, 300}, {300, 500}, {100, 300}
                };

                auto result = ClipPolygonSutherlandHodgman(star, clip);

                std::vector<Point2D> star2d, clip2d, res2d;
                for (auto &p: star) star2d.push_back(p.ToPoint2D());
                for (auto &p: clip) clip2d.push_back(p.ToPoint2D());
                for (auto &p: result) res2d.push_back(p.ToPoint2D());

                std::ignore = DrawPolygonEdges(img, star2d, {100, 100, 100});
                std::ignore = DrawPolygonEdges(img, clip2d, {255, 255, 0});

                if (!res2d.empty()) {
                    FillPolygonEvenOdd(img, res2d, {0, 100, 0});
                    std::ignore = DrawPolygonEdges(img, res2d, {0, 255, 255});
                }
            }, 600, 600);

            res = DoSmthWithImage("lab4_edge_alignment_test.png", [](Image &img) {
                img.ResizeRgb(400, 400);
                std::vector<Point> rect = {{100, 100}, {300, 100}, {300, 300}, {100, 300}};

                Point p0(50, 100), p1(350, 100);

                std::vector<Point2D> rect2d;
                for (auto &p: rect) rect2d.push_back(p.ToPoint2D());
                std::ignore = DrawPolygonEdges(img, rect2d, {255, 255, 0});

                std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {60, 60, 60});

                if (CyrusBeckClipSegment(p0, p1, rect)) {
                    std::ignore = DrawLine(img, p0.ToPoint2D(), p1.ToPoint2D(), {255, 0, 0});
                }
            }, 400, 400);

            return {};
        },
        [](const Laboratory5 &c) -> std::expected<void, boost::system::error_code> {
            Point3 center(0, 0, 0);
            Point3 size(150, 100, 80);
            Point3 axis(c.axis_x, c.axis_y, c.axis_z);

            for (int i = 0; i < c.frames; ++i) {
                auto angle = (2.0 * M_PI * i) / c.frames;

                auto fname_par = fmt::format("{}_parallel_{:03d}.png", c.output_prefix, i);
                auto res = DoSmthWithImage(fname_par, [&](Image &img) {
                    RenderParallelepiped(img, center, size, axis, angle, ProjectionType::kParallel, c.k);
                }, 400, 400);

                if (!res) {
                    return std::unexpected{res.error()};
                }

                auto fname_persp = fmt::format("{}_perspective_{:03d}.png", c.output_prefix, i);
                res = DoSmthWithImage(fname_persp, [&](Image &img) {
                    RenderParallelepiped(img, center, size, axis, angle, ProjectionType::kPerspective, c.k);
                }, 400, 400);

                if (!res) {
                    return std::unexpected{res.error()};
                }
            }

            for (int i = 0; i < c.frames; ++i) {
                double t = (2.0 * M_PI * i) / c.frames;

                double radius = 150.0;
                Point3 pos1(radius * cos(t), 0, radius * sin(t));
                Point3 pos2(radius * cos(t + M_PI), 0, radius * sin(t + M_PI));

                SceneObject cuboid1{
                    .center = pos1,
                    .size = {60, 60, 60},
                    .rotation_axis = {0, 1, 0},
                    .rotation_angle = t * 2,
                    .face_colors = {{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, {255, 0, 255}, {0, 255, 255}}
                };

                SceneObject cuboid2{
                    .center = pos2,
                    .size = {80, 40, 50},
                    .rotation_axis = {1, 1, 0},
                    .rotation_angle = -t * 1.5,
                    .face_colors = {{128, 0, 0}, {0, 128, 0}, {0, 0, 128}, {128, 128, 0}, {128, 0, 128}, {0, 128, 128}}
                };

                auto fname = fmt::format("{}_{:03d}.png", c.output_prefix, i);
                auto res = DoSmthWithImage(fname, [&](Image &img) {
                    RenderLab5Scene(img,
                                    cuboid1, PerspectiveType::kTwoPoint,
                                    cuboid2, PerspectiveType::kThreePoint,
                                    c.k);
                }, 600, 600);

                if (!res) return std::unexpected{res.error()};
            }

            return {};
        },
        [](const Homework1 &c) -> std::expected<void, boost::system::error_code> {
            std::vector<Point2D> poly = {
                {100, 100}, {400, 100}, {100, 400}, {400, 400}
            };

            auto res = DoSmthWithImage(c.output_filename + "_original.png", [&](Image &img) {
                DrawPolygonEdges(img, poly, {100, 100, 100});
            }, 500, 500);

            if (!res) {
                return std::unexpected{res.error()};
            }

            auto contour = ExtractExteriorContour(poly, 500, 500);

            res = DoSmthWithImage(c.output_filename, [&](Image &img) {
                FillPolygonNonZero(img, poly, {50, 50, 150});
                DrawPolygonEdges(img, contour, {255, 255, 0});
            }, 500, 500);

            if (!res) {
                return std::unexpected{res.error()};
            }

            fmt::print("Extracted contour size: {} points\n", contour.size());

            res = DoSmthWithImage("hw1_contour.png", [](Image &img) {
                std::vector<Point2D> star = {{50, 50}, {450, 450}, {450, 50}, {50, 450}, {250, 10}};

                FillPolygonNonZero(img, star, {50, 0, 0});
                auto contour = ExtractExteriorContour(star, 500, 500);

                for (size_t i = 0; i < contour.size(); ++i) {
                    Point2D p1 = contour[i];
                    Point2D p2 = contour[(i + 1) % contour.size()];

                    DrawLine(img, p1, p2, {0, 255, 0});
                }
            }, 500, 500);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("hw1_complex_star.png", [](Image &img) {
                img.ResizeRgb(600, 600);
                std::vector<Point2D> star = {
                    {300, 50}, {350, 550}, {100, 200}, {500, 200},
                    {250, 550}, {300, 50}, {450, 400}, {150, 400}
                };

                FillPolygonNonZero(img, star, {40, 40, 80});

                auto contour = ExtractExteriorContour(star, 600, 600);

                for (size_t i = 0; i < contour.size(); ++i) {
                    DrawThickLine(img, contour[i], contour[(i + 1) % contour.size()], 2.0, {255, 0, 0});
                }
                fmt::print("HW1 Complex Star: Contour size {} points\n", contour.size());
            }, 600, 600);

            res = DoSmthWithImage("hw1_bowtie.png", [](Image &img) {
                img.ResizeRgb(500, 500);
                std::vector<Point2D> poly = {{100, 100}, {400, 400}, {400, 100}, {100, 400}};

                auto contour = ExtractExteriorContour(poly, 500, 500);

                DrawPolygonEdges(img, poly, {80, 80, 80});

                for (size_t i = 0; i < contour.size(); ++i) {
                    DrawThickLine(img, contour[i], contour[(i + 1) % contour.size()], 3.0, {0, 255, 0});
                }
                fmt::print("Test Bowtie: Original 4 points -> Contour {} points\n", contour.size());
            }, 500, 500);

            res = DoSmthWithImage("hw1_star_silhouette.png", [](Image &img) {
                img.ResizeRgb(500, 500);
                std::vector<Point2D> star = {
                    {250, 50}, {310, 450}, {50, 200}, {450, 200}, {190, 450}
                };

                auto contour = ExtractExteriorContour(star, 500, 500);

                FillPolygonNonZero(img, star, {40, 40, 60});
                DrawPolygonEdges(img, star, {100, 100, 100});

                for (size_t i = 0; i < contour.size(); ++i) {
                    DrawLine(img, contour[i], contour[(i + 1) % contour.size()], {255, 255, 0});
                }
                fmt::print("Test Star: Original 5 points -> Contour {} points\n", contour.size());
            }, 500, 500);

            res = DoSmthWithImage("hw1_overlap_mesh.png", [](Image &img) {
                img.ResizeRgb(600, 600);
                std::vector<Point2D> complex = {
                    {100, 100}, {500, 100}, {500, 500}, {200, 500},
                    {200, 50}, {400, 50}, {400, 550}, {50, 550}
                };

                auto contour = ExtractExteriorContour(complex, 600, 600);

                DrawPolygonEdges(img, complex, {70, 70, 70});

                for (size_t i = 0; i < contour.size(); ++i) {
                    DrawThickLine(img, contour[i], contour[(i + 1) % contour.size()], 2.0, {255, 0, 0});
                }
                fmt::print("Test Overlap: Original 8 points -> Contour {} points\n", contour.size());
            }, 600, 600);

            return {};
        },

        [](const Homework2 &c) -> std::expected<void, boost::system::error_code> {
            Point center(c.center_x, c.center_y);
            auto arc_points = MakeCircleArc(center, c.radius, c.angle_start, c.angle_end);

            auto res = DoSmthWithImage(c.output_filename, [&](Image &img) {
                for (const auto &p: arc_points) {
                    auto p2d = p.ToPoint2D();
                    img.GetRGBPixel(p2d.x, p2d.y) = {0, 255, 0};
                }

                DrawLine(img, center.ToPoint2D(), arc_points.front().ToPoint2D(), {100, 100, 100});
                DrawLine(img, center.ToPoint2D(), arc_points.back().ToPoint2D(), {100, 100, 100});
            }, 600, 600);

            if (!res) {
                return std::unexpected{res.error()};
            }

            res = DoSmthWithImage("hw2_arc_gallery.png", [](Image &img) {
                auto a1 = MakeCircleArc({150, 150}, 100, 0, 90);
                auto a2 = MakeCircleArc({350, 150}, 100, 0, 360);
                auto a3 = MakeCircleArc({150, 350}, 80, 45, 315);

                for (const auto &p: a1) {
                    auto p2d = p.ToPoint2D();
                    img.GetRGBPixel(p2d.x, p2d.y) = {0, 255, 0};
                }

                for (const auto &p: a2) {
                    auto p2d = p.ToPoint2D();
                    img.GetRGBPixel(p2d.x, p2d.y) = {0, 255, 0};
                }

                for (const auto &p: a3) {
                    auto p2d = p.ToPoint2D();
                    img.GetRGBPixel(p2d.x, p2d.y) = {0, 255, 0};
                }
            }, 600, 600);

            if (!res) {
                return std::unexpected{res.error()};
            }

            return {};
        },

        [](const Homework3 &c) -> std::expected<void, boost::system::error_code> {
            auto image_res = ReadImage(c.input_filename);
            if (!image_res) return std::unexpected{image_res.error()};

            fmt::print("Starting K-Means quantization (k={})...\n", c.k_colors);
            ColorQuantizationKMeans(image_res.value(), c.k_colors);

            auto save = SaveImage(c.output_filename, image_res.value());
            if (!save) return std::unexpected{save.error()};

            fmt::print("Saved quantized image: {}\n", c.output_filename);
            return {};
        }
    };

    try {
        if (cmd.has_value()) {
            auto result = std::visit(v, cmd.value());

            if (!result) {
                fmt::print("Failed with error: {}\n", result.error().message());
                return 1;
            }
        } else {
            std::cout << parser << std::endl;
        }
    } catch (const std::exception &ex) {
        std::cerr << "Fatal: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
