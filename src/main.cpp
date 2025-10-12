#include <algorithm>
#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <memory>
#include <atomic>
#include <mutex>
#include <string>
#include <vector>

#include <boost/system/error_code.hpp>
#include <boost/system/error_code.hpp>
#include <fmt/core.h>
#include <argparse/argparse.hpp>
#include <png.h>

int main(int argc, char *argv[]) {
    std::string outputFile = argc > 1 ? argv[1] : "processes.json";
    try {
        boost::system::error_code ec;
        std::cout << "Boost.System error code: " << ec.message() << "\n";

        fmt::print("Hello from fmt!\n");

        argparse::ArgumentParser program("test");
        program.add_argument("--value").help("Just a test");

        std::cout << "Libpng version: " << PNG_LIBPNG_VER_STRING << "\n";
    } catch (const std::exception &ex) {
        std::cerr << "Fatal error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
