#include <expected>
#include <commands.hpp>
#include <lab1.hpp>
#include <lab2.hpp>
#include <image.hpp>

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
    command_lab1.add_description("starts laboratory 1");

    argparse::ArgumentParser command_lab1_circle_gray("circle-gray");
    command_lab1_circle_gray.add_description("prints circle gray");
    AddOutputArgument(command_lab1_circle_gray, output);

    argparse::ArgumentParser command_lab1_blending("blend");
    command_lab1_blending.add_description("blending two images");
    AddOutputArgument(command_lab1_blending, output);
    AddInputArgument(command_lab1_blending, input_first, "--input-first");
    AddInputArgument(command_lab1_blending, input_second, "--input-second");
    AddInputArgument(command_lab1_blending, input_third, "--input-alpha");

    command_lab1.add_subparser(command_lab1_circle_gray);
    command_lab1.add_subparser(command_lab1_blending);
    parser.add_subparser(command_lab1);

    argparse::ArgumentParser command_lab2("lab2");
    command_lab2.add_description("starts laboratory 2");
    AddInputArgument(command_lab2, input_first);
    AddOutputArgument(command_lab2, output);

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
            return Laboratory2{};
        }

        return std::nullopt;
    };

    std::optional<Command> cmd = initialize_cmd();

    const auto v = overloaded{
        [](const Laboratory1GrayCircle &c) -> std::expected<void, boost::system::error_code> {
            auto image = Image{c.width, c.height};
            MakeCircularGrayscale(image, c.radius_fraction);
            SaveGrayPng(c.output_filename, image);

            fmt::print("Saved grayscale png image: {}, with width {} and height {}, with radius fraction: {}\n",
                       c.output_filename,
                       c.width, c.height, c.radius_fraction);

            return {};
        },
        [](const Laboratory1Blend &c) -> std::expected<void, boost::system::error_code> {
            const Image A = ReadGrayPng(c.first_filename);
            const Image B = ReadGrayPng(c.second_filename);
            const Image Alpha = ReadGrayPng(c.alpha_filename);

            Image blended = Blend(A, B, Alpha);
            WriteGrayPng(c.output_filename, blended);

            fmt::print("Saved blended png image: {}\n",
                       c.output_filename);

            return {};
        },
        [](const Laboratory2 &c) -> std::expected<void, boost::system::error_code> {
            fmt::print("Unimplemented!\n");
            return std::unexpected{
                boost::system::errc::make_error_code(
                    boost::system::errc::not_supported
                )
            };
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
            std::cout << command_lab1 << std::endl;
            std::cout << command_lab2 << std::endl;
        }
    } catch (const std::exception &ex) {
        std::cerr << "Fatal error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
