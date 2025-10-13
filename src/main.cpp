#include <expected>
#include <commands.hpp>
#include <improcessing.hpp>
#include <image.hpp>
#include <iostream>
#include <filesystem>

#include <boost/system/error_code.hpp>
#include <fmt/core.h>
#include <argparse/argparse.hpp>

auto AddInputArgument(argparse::ArgumentParser &parser, std::string &input,
                      const std::string &name = "--input") -> void {
    parser.add_argument(name)
            .required()
            .store_into(input)
            .help("specify the input file");
}

auto AddOutputArgument(argparse::ArgumentParser &parser, std::string &output) -> void {
    parser.add_argument("-o", "--output")
            .required()
            .store_into(output)
            .help("specify the output file");
}

int main(int argc, char *argv[]) {
    std::string input_first{};
    std::string input_second{};
    std::string input_third{};
    std::string output{};

    argparse::ArgumentParser parser("image-processing", "0.1.0", argparse::default_arguments::none);
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
    AddInputArgument(command_lab2, input_first);
    AddOutputArgument(command_lab2, output);
    command_lab2.add_argument("-n", "--n-levels").required().store_into(n_levels).help("set level for scattering");

    parser.add_subparser(command_lab2);

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
                    .width = 512,
                    .height = 512,
                    .radius_fraction = 0.45,
                    .output_filename = output
                };
            } else if (command_lab1.is_subcommand_used(command_lab1_blending)) {
                return Laboratory1Blend{
                    .first_filename = input_first,
                    .second_filename = input_second,
                    .alpha_filename = input_third,
                    .output_filename = output
                };
            }
        } else if (parser.is_subcommand_used(command_lab2)) {
            return Laboratory2{
                .input_filename = input_first,
                .output_filename = output,
                .n_levels = n_levels
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

            fmt::print("Saved grayscale png image: {}, with width {} and height {}, with radius fraction: {}\n",
                       (std::filesystem::current_path() / c.output_filename).string(),
                       c.width, c.height, c.radius_fraction);

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

            fmt::print("Saved blended png image: {}\n",
                       (std::filesystem::current_path() / c.output_filename).string());

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

            auto save = SaveImage(c.output_filename, result.value());
            if (!save) {
                return std::unexpected{save.error()};
            }

            fmt::print("Saved scatter image: {}, with n_levels: {}\n",
                       (std::filesystem::current_path() / c.output_filename).string(),
                       c.n_levels);

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
