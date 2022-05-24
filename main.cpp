//
// Created by fbasatemur on 4/23/2022.
//

#include "Config.h"
#include "Image.h"
#include "Container.h"
#include "NSST.h"

int main()
{
    long size, bmp_offset;
    int width, height;
    BYTE* image = nullptr;
    NSST nsst;

    // LoadBMP can read only 24 bit image depth
    image = LoadBMP(width, height, size, static_cast<LPCTSTR>(Conf::image_path));

    // Low and High Frequences Coefficients -- 2D Laplacian Pyramid filters
    nsst.AtrousFilters();      // filters is in NSST

    // New filter coefficients are obtained	-- Only once -- Optional
    //filters->mats[1] = Conv2(filters->mats[1], filters->mats[0], "same");
    //filters->mats[3] = Conv2(filters->mats[3], filters->mats[2], "same");

    nsst.ShearingFiltersMyer();

    if (image)
    {
        // Intensity image form
        float* intensity = ConvertBMPToIntensity(image, width, height);

        Tensor i_image;
        i_image.Set(height, width);
        i_image._mat = Eigen::Map<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(intensity, height, width);                        // I => Intensity


        // NSST - Non Subsampled Shearlet Transform
        Cont* dst = nsst.NsstDec1e(&i_image);
        //	INFO
        //  dst->_conts[_cont_num][_cont_depth]
        //  dst->_conts[0][0]			=> AFK is 1 piece and deep	 => 1
        //  dst->_conts[1..4][deep]	    => YFK is 4 pieces and deeps => {8, 8, 16, 16}


        // Inverse NSST
        Tensor* inverse = nsst.NsstRec1(dst);


        // save Inverse NSST result
        BYTE* inverse_bmp = ConvertIntensityToBMP(inverse->_mat, bmp_offset);
        SaveBMP(inverse_bmp, inverse->_w, inverse->_h,bmp_offset, static_cast<LPCTSTR>(Conf::save_image_path));

        delete inverse_bmp;
        delete dst;
        delete inverse;
        delete intensity;
        delete image;
    }

    return 0;
}

