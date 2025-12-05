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
                // {11, 25}, {30, 40}, {26, 20}, {18, 16}, {12, 26}, {10, 14}
                {0, 0}, {0, 20}, {20, 20}, {20, 0}
                // {50, 0},
                // {60, 35},
                // {100, 35},
                // {65, 55},
                // {80, 100},
                // {50, 70},
                // {20, 100},
                // {35, 55},
                // {0, 35},
                // {40, 35}
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
