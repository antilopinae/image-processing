#pragma once

#include <string>

class Image;

auto SaveGrayPng(const std::string &filename, const Image &image) -> void;

auto MakeCircularGrayscale(const Image &image, double radius_fraction = 0.45) -> void;

auto ReadGrayPng(const std::string& filename) -> Image;

auto WriteGrayPng(const std::string& filename, const Image& img) -> void;

auto Blend(const Image& A, const Image& B, const Image& Alpha)-> Image;

