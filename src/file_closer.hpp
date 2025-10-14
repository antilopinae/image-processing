#pragma once

#include <cstdio>

namespace improcessing
{
    namespace details
    {
        struct FileCloser
        {
            void operator()(FILE* f) const noexcept
            {
                if (f)
                    std::fclose(f);
            }
        };
    } // namespace details
} // namespace improcessing
