#pragma once

#include <cstdint>
#include <string>
#include <variant>

template<class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

struct Laboratory1GrayCircle {
    const std::uint64_t width, height;
    const double radius_fraction;
    const std::string &output_filename;
};

struct Laboratory1Blend {
    const std::string &first_filename;
    const std::string &second_filename;
    const std::string &alpha_filename;
    const std::string &output_filename;
};

struct Laboratory2 {
    const std::string &input_filename;
    const std::string &output_filename;
    const uint8_t n_levels;
};

using Command = std::variant<Laboratory1GrayCircle, Laboratory1Blend, Laboratory2>;
