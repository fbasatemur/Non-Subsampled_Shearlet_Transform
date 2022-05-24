#ifndef NSST_IMAGE_H
#define NSST_IMAGE_H

#include <Eigen/Core>
#include <Windows.h>

BYTE* LoadBMP(int& width, int& height, long& size, LPCTSTR bmpfile);
float* ConvertBMPToIntensity(BYTE* Buffer, int width, int height);
BYTE* ConvertIntensityToBMP(const Eigen::ArrayXXf& Buffer, long& newsize);
bool SaveBMP(BYTE* Buffer, int width, int height, long paddedsize, LPCTSTR bmpfile);
float* ConvertBMPToYCbCr(BYTE* Buffer, int width, int height);
#endif