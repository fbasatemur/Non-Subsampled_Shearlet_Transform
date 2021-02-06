#pragma once
#include <Windows.h>
struct imagener {
	float* Real;
	float* Im;
	int w;
	int h;
};
void FFT2D(float* img, float* out_real, float* out_imag, int width, int height);
void IFFT2D(float* Output_real, float* Output_imag, float* Input_real, float* Input_imag, int width, int height);