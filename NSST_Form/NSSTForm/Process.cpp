// 2-D FFT kodu (Tek boyutlu FFT yardimiyla olusturuldu).
#include "Process.h"
#include "FFT.h"

void FFT2D(float* img, float* Output_real, float* Output_img, int width, int height)
{
	int i, j;
	float* input_real = new float[width];
	float* input_im = new float[width];
	float* out_real = new float[width];
	float* out_im = new float[width];
	float* Real = new float[width * height];
	float* Im = new float[width * height];

	for (j = 0; j < width; j++)
		input_im[j] = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++)
			input_real[j] = img[i * width + j];

		fft(width, input_real, input_im, out_real, out_im);
		for (j = 0; j < width; j++) {
			Real[i * width + j] = out_real[j];
			Im[i * width + j] = out_im[j];
		}
	}
	delete[] input_im;
	delete[] input_real;
	delete[] out_real;
	delete[] out_im;
	input_real = new float[height];
	input_im = new float[height];
	out_real = new float[height];
	out_im = new float[height];
	for (j = 0; j < width; j++) {
		for (i = 0; i < height; i++) {
			input_real[i] = Real[i * width + j];
			input_im[i] = Im[i * width + j];
		}

		fft(height, input_real, input_im, out_real, out_im);
		for (i = 0; i < height; i++) {
			Output_real[i * width + j] = out_real[i];
			Output_img[i * width + j] = out_im[i];
		}
	}
	delete[] input_im;
	delete[] input_real;
	delete[] out_real;
	delete[] out_im;
	delete[] Real;
	delete[] Im;
}//FFT2D

// Inverse FFT kodu
void IFFT2D(float* Output_real, float* Output_imag, float* Input_real, float* Input_imag, int width, int height)
{
	int i, j;
	float* In_re, * In_im, * O_re, * O_im;

	In_re = new float[width];
	In_im = new float[width];
	O_re = new float[width];
	O_im = new float[width];
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			In_re[j] = Input_real[i * width + j];
			In_im[j] = Input_imag[i * width + j];
		}
		ifft(width, In_re, In_im, O_re, O_im);
		for (j = 0; j < width; j++) {
			Output_real[i * width + j] = O_re[j];
			Output_imag[i * width + j] = O_im[j];
		}
	}
	delete[] In_re; delete[] In_im; delete[] O_re; delete[] O_im;
	In_re = new float[height];
	In_im = new float[height];
	O_re = new float[height];
	O_im = new float[height];

	for (j = 0; j < width; j++)
	{
		for (i = 0; i < height; i++)
		{
			In_re[i] = Output_real[i * width + j];
			In_im[i] = Output_imag[i * width + j];
		}
		ifft(height, In_re, In_im, O_re, O_im);
		for (i = 0; i < height; i++)
		{
			Output_real[i * width + j] = O_re[i];
			Output_imag[i * width + j] = O_im[i];
		}
	}
	delete[] In_re; delete[] In_im; delete[] O_re; delete[] O_im;
}//ifft2D


