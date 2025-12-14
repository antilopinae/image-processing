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

    auto k_param = 1000.0;
    auto frames = 36;

    command_lab5.add_argument("-k", "--perspective-k").scan<'g', double>().store_into(k_param).default_value(500.0).
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
