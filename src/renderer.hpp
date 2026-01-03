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
    Matrix4x4 CalculateProjectionMatrix(const SceneObject& obj, const Camera& cam)
    {
        auto P = Matrix4x4::Identity();

        auto r = (cam.focal_length > 0) ? (1.0 / cam.focal_length) : 0.0;

        if (obj.type == PerspectiveType::kOnePoint) {
            P.m[2][3] = r;
        } else if (obj.type == PerspectiveType::kTwoPoint) {
            P.m[0][3] = r;
            P.m[1][3] = r;
        } else if (obj.type == PerspectiveType::kThreePoint) {
            P.m[0][3] = r;
            P.m[1][3] = r;
            P.m[2][3] = r;
        }
        return P;
    }

    Matrix4x4 CalculateModelMatrix(const SceneObject& obj)
    {
        auto T = Matrix4x4::Translation(obj.center.x, obj.center.y, obj.center.z);

        const auto rx = Matrix4x4::RotationX(obj.rotation.x);
        const auto ry = Matrix4x4::RotationY(obj.rotation.y);
        const auto rz = Matrix4x4::RotationZ(obj.rotation.z);

        auto R = rz * ry * rx;

        return T * R;
    }

    void Render(Image& img,
                const Cube& cube,
                const SceneObject& settings,
                const Camera& cam,
                const std::vector<Light>& lights,
                std::vector<double>& z_buffer)
    {
        auto modelMatrix = CalculateModelMatrix(settings);
        auto viewM       = Matrix4x4::Translation(-cam.position.x, -cam.position.y, -cam.position.z);
        auto projM       = CalculateProjectionMatrix(settings, cam);

        for (const auto& face : cube.faces) {
            std::vector<VertexData> faceVertices;

            for (auto i = 0; i < 4; ++i) {
                VertexData v;
                auto local_p = cube.vertices[face.indices[i]];

                v.world_pos = modelMatrix.Transform(local_p);
                faceVertices.push_back(v);
            }

            auto edge1       = faceVertices[1].world_pos - faceVertices[0].world_pos;
            auto edge2       = faceVertices[2].world_pos - faceVertices[0].world_pos;
            auto face_normal = edge1.Cross(edge2).Normalized();

            for (auto& v : faceVertices) {
                v.normal = face_normal;
            }

            auto clipped = ClipPolygonToNearPlane(faceVertices, cam.position.z + cam.near_plane);
            if (clipped.size() < 3) {
                continue;
            }

            for (auto& v : clipped) {
                auto p_view = v.world_pos - cam.position;

                auto X = p_view.x;
                auto Y = p_view.y;
                auto Z = p_view.z;
                auto H = X * projM.m[0][3] + Y * projM.m[1][3] + Z * projM.m[2][3] + 1.0;

                // x' = X/H, y' = Y/H
                auto x_norm = X / H;
                auto y_norm = Y / H;

                v.screen_pos.x = (img.WidthRgb() / 2.0) + x_norm * cam.focal_length;
                v.screen_pos.y = (img.HeightRgb() / 2.0) - y_norm * cam.focal_length;
                v.screen_pos.z = Z;
            }

            for (size_t i = 1; i < clipped.size() - 1; ++i) {
                RasterizeTriangle(img,
                                  z_buffer,
                                  clipped[0],
                                  clipped[i],
                                  clipped[i + 1],
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

                        auto p = v0.world_pos * bary.x + v1.world_pos * bary.y + v2.world_pos * bary.z;

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
};

} // namespace improcessing