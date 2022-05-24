//
// Created by fbasatemur on 4/23/2022.
//

#ifndef NSST_DEFINITIONS_H
#define NSST_DEFINITIONS_H

#if __unix__
    #define BYTE uint8_t
    #define UINT uint32_t
#else
    typedef unsigned char BYTE;
    typedef unsigned int UINT;
#endif

typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXf;

#endif //NSST_DEFINITIONS_H
