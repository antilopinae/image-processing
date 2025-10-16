#pragma once

namespace improcessing {
    template<typename T>
    struct PixelRGB {
        union {
            struct {
                T r;
                T g;
                T b;
            };

            T data[3];
        };

        explicit PixelRGB() : r(0), g(0), b(0) {
        }

        PixelRGB(T r_, T g_, T b_) : r(r_), g(g_), b(b_) {
        }
    };
}
