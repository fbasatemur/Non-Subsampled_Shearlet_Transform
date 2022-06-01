//
// Created by fbasatemur on 4/23/2022.
//

#include "Config.h"
#include "Image.h"
#include "Container.h"
#include "NSST.h"

int main()
{
    ImgBMP img;
    NSST nsst;

    // LoadBMP can read only 24 bit image depth
    img.LoadImage(Conf::image_path);
    int width = (int)img.info_header_data.biWidth, height = (int)img.info_header_data.biHeight;

    // Low and High Frequences Coefficients -- 2D Laplacian Pyramid filters
    nsst.AtrousFilters();      // filters is in NSST

    // New filter coefficients are obtained	-- Only once -- Optional
    //filters->mats[1] = Conv2(filters->mats[1], filters->mats[0], "same");
    //filters->mats[3] = Conv2(filters->mats[3], filters->mats[2], "same");

    nsst.ShearingFiltersMyer();

    if (img.bmp_p)
    {
        // Intensity image form
        img.BMP2Intensity();

        Tensor i_image;
        i_image.Set(height, width);
        i_image._mat = Eigen::Map<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(img.intensity_p, height, width);


        // NSST - Non Subsampled Shearlet Transform
        Cont* dst = nsst.Dec(&i_image);
        //	INFO
        //  dst->_conts[_cont_num][_cont_depth]
        //  dst->_conts[0][0]			=> AFK is 1 piece and deep	 => 1
        //  dst->_conts[1..4][deep]	    => YFK is 4 pieces and deeps => {8, 8, 16, 16}


        // Inverse NSST
        Tensor* inverse = nsst.Rec(dst);


        // Inverse NSST result are saved
        img.Intensity2BMP(inverse->_mat.data());
        img.SaveImage(Conf::save_image_path);

        delete dst;
        delete inverse;
    }

    return 0;
}

