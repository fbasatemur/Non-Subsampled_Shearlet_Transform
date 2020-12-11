// 2-D FFT kodu (Tek boyutlu FFT yardimiyla olusturuldu).
#include "Process.h"
#include "FFT.h"

void FFT2D(double* img, double* Output_real, double* Output_img, int width, int height)
{
	int i, j;
	double* input_real = new double[width];
	double* input_im = new double[width];
	double* out_real = new double[width];
	double* out_im = new double[width];
	double* Real = new double[width * height];
	double* Im = new double[width * height];

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
	input_real = new double[height];
	input_im = new double[height];
	out_real = new double[height];
	out_im = new double[height];
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
void IFFT2D(double* Output_real, double* Output_imag, double* Input_real, double* Input_imag, int width, int height)
{
	int i, j;
	double* In_re, * In_im, * O_re, * O_im;

	In_re = new double[width];
	In_im = new double[width];
	O_re = new double[width];
	O_im = new double[width];
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
	In_re = new double[height];
	In_im = new double[height];
	O_re = new double[height];
	O_im = new double[height];

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


