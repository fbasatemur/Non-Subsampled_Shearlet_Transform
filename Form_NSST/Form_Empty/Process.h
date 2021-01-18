#pragma once
#include <Windows.h>
struct imagener {
	double* Real;
	double* Im;
	int w;
	int h;
};
void FFT2D(double* img, double* out_real, double* out_imag, int width, int height);
void IFFT2D(double* Output_real, double* Output_imag, double* Input_real, double* Input_imag, int width, int height);