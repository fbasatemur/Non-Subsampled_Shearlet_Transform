#ifndef NSST_IMAGE_H
#define NSST_IMAGE_H

#include <Eigen/Core>
#include <vector>
#include <stdio.h>
#include <stdint.h>

typedef struct
{
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} __attribute__((__packed__)) BitMapFileHeader;

typedef struct
{
    uint32_t biSize;
    int32_t biWidth;
    int32_t biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t biXPelsPerMeter;
    int32_t biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} __attribute__((__packed__)) BitMapInfoHeader;

typedef struct
{
    uint8_t blue;
    uint8_t green;
    uint8_t red;
    uint8_t reserved;
} __attribute__((__packed__)) BGRMap;

struct ImgBMP
{
public:
    BitMapFileHeader file_header_data;
    BitMapInfoHeader info_header_data;
    BGRMap color_map[256];
    uint8_t* bmp_p = nullptr;
    float* intensity_p = nullptr;

    bool LoadImage(const char* path);
    void BMP2Intensity();
    void Intensity2BMP(const float* intensity);
    bool SaveImage(const char* path);

    ~ImgBMP(){
        delete bmp_p;
        delete intensity_p;
    }
}__attribute__((__packed__));

bool ImgBMP::LoadImage(const char* path)
{
    FILE* pFile = fopen(path, "rb");
    if (!pFile)
    {
        return 0;
    }
    // Processing
    fread(&file_header_data, sizeof(BitMapFileHeader), 1, pFile);
    if (file_header_data.bfType == 0x4D42) // Check is it an RGB file
    {
        // Get Channel num of a pixel
        int channels = 0;
        fread(&info_header_data, sizeof(BitMapInfoHeader), 1, pFile);
        if (info_header_data.biBitCount == 8)// grayscale format
        {

            channels = 1;
            fread(&color_map, sizeof(BGRMap), 256, pFile);
        }
        else if (info_header_data.biBitCount == 24)// RGB format
        {
            channels = 3;
        }

        // Get offset of every scanline,length(scanline)=length(pixel)+offset
        int offset = 0;
        int linelength = info_header_data.biWidth * channels;
        offset = linelength % 4;
        if (offset > 0)
        {
            offset = 4 - offset;
        }

        // Read Pixel
        bmp_p = (uint8_t*)malloc(info_header_data.biHeight * linelength * sizeof(uint8_t));
        for (int i = 0; i < info_header_data.biHeight; i++)
        {
            fread(bmp_p + i * linelength, linelength, 1, pFile);
            fseek(pFile, offset, SEEK_CUR);
        }
    }
    else
    {
        return false;
    }

    fclose(pFile);
    return true;
}

void ImgBMP::BMP2Intensity(){

    intensity_p = new float[info_header_data.biWidth * info_header_data.biHeight];
    long bufpos = 0;
    for (int i = 0; i < info_header_data.biWidth * info_header_data.biHeight; ++i, bufpos+=3){
        intensity_p[i] = (float)(int)(0.11 * bmp_p[bufpos + 2] + 0.59 * bmp_p[bufpos + 1] + 0.3 * bmp_p[bufpos]);

        /*
            Buffer[bufpos]	   => Blue
            Buffer[bufpos + 1] => Green
            Buffer[bufpos + 2] => Red
        */
    }
}

void ImgBMP::Intensity2BMP(const float* intensity) {

    long bufpos = 0;
    for (int i = 0; i < info_header_data.biWidth * info_header_data.biHeight; ++i, bufpos+=3){
            bmp_p[bufpos] = bmp_p[bufpos + 1] = bmp_p[bufpos + 2] = intensity[i];
    }
}

bool ImgBMP::SaveImage(const char* path)
{
    FILE* pFile = fopen(path, "wb");
    if (!pFile)
    {
        return 0;
    }

    // Processing
    fwrite(&file_header_data, sizeof(BitMapFileHeader), 1, pFile);
    fwrite(&info_header_data, sizeof(BitMapInfoHeader), 1, pFile);
    // Get Channel num of a pixel
    int channels = 0;
    if (info_header_data.biBitCount == 8)
    {
        channels = 1;
        fwrite(&color_map, sizeof(BGRMap), 256, pFile);
    }
    else if (info_header_data.biBitCount == 24)
    {
        channels = 3;
    }
    // Get offset of every scanline,length(scanline)=length(pixel)+offset
    int offset = 0;
    int linelength = info_header_data.biWidth * channels;
    offset = (channels * info_header_data.biWidth) % 4;
    if (offset > 0)
    {
        offset = 4 - offset;
    }
    // Write Pixel
    uint8_t pixVal = 0;
    for (int i = 0; i < info_header_data.biHeight; i++)
    {
        fwrite(bmp_p + i * linelength, sizeof(uint8_t), linelength, pFile);
        fwrite(&pixVal, sizeof(uint8_t), offset, pFile);
    }
    fclose(pFile);
    return true;
}

#endif