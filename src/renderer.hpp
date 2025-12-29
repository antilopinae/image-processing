#pragma once

#include <camera.hpp>
#include <cube.hpp>
#include <image.hpp>
#include <improcessing.hpp>
#include <light.hpp>
#include <material.hpp>
#include <matrix4x4.hpp>

#include <scene_object.hpp>

namespace improcessing {

class Renderer {
public:
    const int INSIDE = 0 << 0; // 0000
    const int LEFT   = 1 << 1; // 0001
    const int RIGHT  = 1 << 2; // 0010
    const int BOTTOM = 1 << 3; // 0100
    const int TOP    = 1 << 4; // 1000

    int ComputeCode(double x, double y, double width, double height)
    {
        int code = INSIDE;
        if (x < 0) {
            code |= LEFT;
        } else if (x >= width) {
            code |= RIGHT;
        }
        if (y < 0) {
            code |= BOTTOM;
        } else if (y >= height) {
            code |= TOP;
        }
        return code;
    }

    void ClippedDrawLine(Image& img, Point p1_raw, Point p2_raw, Pixel color)
    {
        double x1 = p1_raw.x, y1 = p1_raw.y;
        double x2 = p2_raw.x, y2 = p2_raw.y;
        double w = (double)img.WidthRgb();
        double h = (double)img.HeightRgb();

        while (true) {
            int code1 = ComputeCode(x1, y1, w, h);
            int code2 = ComputeCode(x2, y2, w, h);

            if ((code1 | code2) == 0) {
                DrawLine(img, {(size_t)x1, (size_t)y1}, {(size_t)x2, (size_t)y2}, color);
                break;
            } else if (code1 & code2) {
                break;
            } else {
                int codeOut = code1 ? code1 : code2;
                double x, y;

                if (codeOut & TOP) {
                    x = x1 + (x2 - x1) * (h - 1 - y1) / (y2 - y1);
                    y = h - 1;
                } else if (codeOut & BOTTOM) {
                    x = x1 + (x2 - x1) * (0 - y1) / (y2 - y1);
                    y = 0;
                } else if (codeOut & RIGHT) {
                    y = y1 + (y2 - y1) * (w - 1 - x1) / (x2 - x1);
                    x = w - 1;
                } else { // LEFT
                    y = y1 + (y2 - y1) * (0 - x1) / (x2 - x1);
                    x = 0;
                }

                if (codeOut == code1) {
                    x1 = x;
                    y1 = y;
                } else {
                    x2 = x;
                    y2 = y;
                }
            }
        }
    }

    void Render(Image& img, const Cube& cube, const SceneObject& settings, const Camera& cam)
    {
        Matrix4x4 modelMatrix = CalculateModelMatrix(settings);

        std::vector<Point> projected_points;

        for (const auto& v : cube.vertices) {
            Point3 p = modelMatrix.Transform(v);
            p.x += settings.center.x;
            p.y += settings.center.y;
            p.z += settings.center.z;

            double z_eff = p.z - cam.position.z;
            if (z_eff < 0.1) {
                z_eff = 0.1;
            }

            double x_scr = (img.WidthRgb() / 2.0) + (p.x * cam.focal_length) / z_eff;
            double y_scr = (img.HeightRgb() / 2.0) - (p.y * cam.focal_length) / z_eff;

            projected_points.push_back({x_scr, y_scr});
        }

        for (const auto& edge : cube.edges) {
            ClippedDrawLine(img, projected_points[edge.start], projected_points[edge.end], settings.color);
        }
    }

    void Render(Image& img,
                const Cube& cube,
                const SceneObject& settings,
                const Camera& cam,
                const std::vector<Light>& lights)
    {
    }

    Matrix4x4 CalculateModelMatrix(const SceneObject& obj)
    {
        Matrix4x4 rot = Matrix4x4::Identity();

        if (obj.type == PerspectiveType::kTwoPoint) {
            rot = Matrix4x4::RotationY(obj.rotation.y);
        } else if (obj.type == PerspectiveType::kThreePoint) {
            const Matrix4x4 rx = Matrix4x4::RotationX(obj.rotation.x);
            const Matrix4x4 ry = Matrix4x4::RotationY(obj.rotation.y);
            const Matrix4x4 rz = Matrix4x4::RotationZ(obj.rotation.z);

            rot = rz * ry * rx;
        }

        return rot;
    }
};

} // namespace improcessing
