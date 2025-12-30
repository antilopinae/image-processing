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

struct VertexData {
    Point3 world_pos;
    Point3 normal;
    Point3 screen_pos;
};

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
                const std::vector<Light>& lights,
                std::vector<double>& z_buffer)
    {
        auto modelMatrix = CalculateModelMatrix(settings);

        for (const auto& face : cube.faces) {
            std::vector<VertexData> faceVertices;

            for (auto i = 0; i < 4; ++i) {
                VertexData v;
                auto local_p = cube.vertices[face.indices[i]];
                v.world_pos  = modelMatrix.Transform(local_p);
                v.world_pos  = v.world_pos + settings.center;

                faceVertices.push_back(v);
            }

            auto edge1       = faceVertices[1].world_pos - faceVertices[0].world_pos;
            auto edge2       = faceVertices[2].world_pos - faceVertices[0].world_pos;
            auto face_normal = edge1.Cross(edge2).Normalized();

            for (auto& v : faceVertices) {
                v.normal = face_normal;
            }

            std::vector<VertexData> clippedPolygon = ClipPolygonToNearPlane(faceVertices,
                                                                            cam.position.z + cam.near_plane);

            if (clippedPolygon.size() < 3) {
                continue;
            }

            for (auto& v : clippedPolygon) {
                auto z_eff     = v.world_pos.z - cam.position.z;
                v.screen_pos.x = (img.WidthRgb() / 2.0) + (v.world_pos.x * cam.focal_length) / z_eff;
                v.screen_pos.y = (img.HeightRgb() / 2.0) - (v.world_pos.y * cam.focal_length) / z_eff;
                v.screen_pos.z = z_eff;
            }

            for (size_t i = 1; i < clippedPolygon.size() - 1; ++i) {
                RasterizeTriangle(img,
                                  z_buffer,
                                  clippedPolygon[0],
                                  clippedPolygon[i],
                                  clippedPolygon[i + 1],
                                  lights,
                                  settings,
                                  cam.position);
            }
        }
    }

    void RasterizeTriangle(Image& img,
                           std::vector<double>& z_buffer,
                           const VertexData& v0,
                           const VertexData& v1,
                           const VertexData& v2,
                           const std::vector<Light>& lights,
                           const SceneObject& obj,
                           const Point3& viewPos)
    {
        auto minX = std::max(0, (int)std::floor(std::min({v0.screen_pos.x, v1.screen_pos.x, v2.screen_pos.x})));
        auto maxX = std::min((int)img.WidthRgb() - 1,
                             (int)std::ceil(std::max({v0.screen_pos.x, v1.screen_pos.x, v2.screen_pos.x})));
        auto minY = std::max(0, (int)std::floor(std::min({v0.screen_pos.y, v1.screen_pos.y, v2.screen_pos.y})));
        auto maxY = std::min((int)img.HeightRgb() - 1,
                             (int)std::ceil(std::max({v0.screen_pos.y, v1.screen_pos.y, v2.screen_pos.y})));

        auto invZ0 = 1.0 / v0.screen_pos.z;
        auto invZ1 = 1.0 / v1.screen_pos.z;
        auto invZ2 = 1.0 / v2.screen_pos.z;

        auto pZ0 = v0.world_pos * invZ0;
        auto pZ1 = v1.world_pos * invZ1;
        auto pZ2 = v2.world_pos * invZ2;

        for (auto y = minY; y <= maxY; ++y) {
            for (auto x = minX; x <= maxX; ++x) {
                auto bary =
                    CalculateBarycentric(Point3{x + 0.5, y + 0.5, 0}, v0.screen_pos, v1.screen_pos, v2.screen_pos);

                if (bary.x >= -1e-7 && bary.y >= -1e-7 && bary.z >= -1e-7) {
                    auto interpolatedInvZ = invZ0 * bary.x + invZ1 * bary.y + invZ2 * bary.z;
                    auto z                = 1.0 / interpolatedInvZ;

                    auto idx = y * img.WidthRgb() + x;
                    if (z < z_buffer[idx]) {
                        z_buffer[idx] = z;

                        auto p = (pZ0 * bary.x + pZ1 * bary.y + pZ2 * bary.z) * z;

                        auto pix              = ComputePhong(p, v0.normal, lights, obj, viewPos);
                        img.GetRGBPixel(x, y) = pix;
                    }
                }
            }
        }
    }

    std::vector<VertexData> ClipPolygonToNearPlane(const std::vector<VertexData>& vertices, double nearZ)
    {
        std::vector<VertexData> result;

        for (auto i = 0u; i < vertices.size(); ++i) {
            const auto& v1 = vertices[i];
            const auto& v2 = vertices[(i + 1) % vertices.size()];

            auto v1_inside = v1.world_pos.z >= nearZ;
            auto v2_inside = v2.world_pos.z >= nearZ;

            if (v1_inside && v2_inside) {
                result.push_back(v2);
            } else if (v1_inside && !v2_inside) {
                auto t = (nearZ - v1.world_pos.z) / (v2.world_pos.z - v1.world_pos.z);
                result.push_back(InterpolateVertex(v1, v2, t));
            } else if (!v1_inside && v2_inside) {
                auto t = (nearZ - v1.world_pos.z) / (v2.world_pos.z - v1.world_pos.z);
                result.push_back(InterpolateVertex(v1, v2, t));
                result.push_back(v2);
            }
        }
        return result;
    }

    VertexData InterpolateVertex(const VertexData& v1, const VertexData& v2, double t)
    {
        VertexData res;
        res.world_pos = v1.world_pos + (v2.world_pos - v1.world_pos) * t;
        res.normal    = (v1.normal + (v2.normal - v1.normal) * t).Normalized();
        return res;
    }

    Pixel ComputePhong(const Point3& p,
                       const Point3& n,
                       const std::vector<Light>& lights,
                       const SceneObject& obj,
                       const Point3& viewPos)
    {
        auto surfaceColor = obj.color.ToPoint3();
        const auto& mat   = obj.material;

        Point3 totalColor(0, 0, 0);
        auto V = (viewPos - p).Normalized();

        for (const auto& light : lights) {
            auto lightColor = light.color.ToPoint3();
            auto L          = (light.position - p).Normalized();

            auto ambient       = lightColor * (light.ambient * mat.k_ambient);
            Point3 ambientPart = {ambient.x * surfaceColor.x, ambient.y * surfaceColor.y, ambient.z * surfaceColor.z};

            auto cosTheta      = std::max(0.0, n.Dot(L));
            auto diffuse       = lightColor * (light.diffuse * mat.k_diffuse * cosTheta);
            Point3 diffusePart = {diffuse.x * surfaceColor.x, diffuse.y * surfaceColor.y, diffuse.z * surfaceColor.z};

            auto R            = (n * (2.0 * n.Dot(L)) - L).Normalized();
            auto cosPhi       = std::max(0.0, R.Dot(V));
            auto specFactor   = std::pow(cosPhi, mat.shininess);
            auto specularPart = lightColor * (light.specular * mat.k_specular * specFactor);

            totalColor = totalColor + ambientPart + diffusePart + specularPart;
        }

        return Pixel::FromPoint3(totalColor);
    }

    Point3 CalculateBarycentric(Point3 p, Point3 a, Point3 b, Point3 c)
    {
        auto v0 = b - a, v1 = c - a, v2 = p - a;
        auto den = v0.x * v1.y - v1.x * v0.y;
        if (std::abs(den) < 1e-8) {
            return {-1, -1, -1};
        }
        auto v = (v2.x * v1.y - v1.x * v2.y) / den;
        auto w = (v0.x * v2.y - v2.x * v0.y) / den;
        auto u = 1.0 - v - w;
        return {u, v, w};
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
