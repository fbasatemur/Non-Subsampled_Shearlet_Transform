//
// Created by fbasatemur on 4/23/2022.
//

#ifndef NEW_NSST_CONFIG_H
#define NEW_NSST_CONFIG_H
#include <array>

// NSST configuration parameters
namespace Conf{

    struct ShearParameters
    {
        static constexpr int dcompSize = 4;        // K => numbers of YFK
        std::array<int, dcompSize> dcomp;
        std::array<int, dcompSize> dsize;
    };

    static char const* image_path = "../images/barbara.bmp";
    static char const* save_image_path = "../images/regenerated.bmp";
    static char const* lpfilt = "maxflat";

    // Shearlet coefficients
    static ShearParameters sp = {.dcomp = {3, 3, 4, 4 },
                                 .dsize = {32, 32, 16, 16 }};

};

#endif //NEW_NSST_CONFIG_H
