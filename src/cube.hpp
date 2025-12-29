#pragma once

#include <vector>

#include <point.hpp>

namespace improcessing {

struct Face {
    int indices[4];
};

class Cube {
public:
    struct Edge {
        int start, end;
    };

    explicit Cube(double size)
    {
        double s = size / 2.0;
        vertices = {
            {-s, -s, -s},
            { s, -s, -s},
            { s,  s, -s},
            {-s,  s, -s},
            {-s, -s,  s},
            { s, -s,  s},
            { s,  s,  s},
            {-s,  s,  s}
        };

        edges = {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 4},
            {0, 4},
            {1, 5},
            {2, 6},
            {3, 7}
        };

        faces = {{{4, 5, 6, 7}}, {{1, 0, 3, 2}}, {{0, 1, 5, 4}}, {{2, 3, 7, 6}}, {{0, 3, 7, 4}}, {{1, 2, 6, 5}}};
    }

    std::vector<Point3> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;
};

} // namespace improcessing