#include <string.h>		//for strcmp
#include "MatlabFuncs.h"
#include <cmath>

# define PI           3.14159265358979323846
#define ERROR (double*)-1

Matrix* Sum(Matrix* mat, int dim) {

	int imageSize = mat->GetSize();

	Matrix* retMat = new Matrix;
	retMat->CreateMatrix(mat->height, mat->width, 1);

	double total = 0.0;

	switch (dim)
	{
	case 1:
		break;
	case 2:
		break;
	case 3:
		for (int i = 0; i < imageSize; i++) {

			total = 0.0;
			for (int d = 0; d < mat->depth; d++)
				total += mat->mat[d * imageSize + i];
			retMat->mat[i] = total;
		}
		
		break;
	default:
		break;
	}

	return retMat;
}

double* Eye(int size) {

	double* identityMat = new double[size * size];

	for (int row = 0; row < size; row++) {
		for (int col = 0; col < size; col++) {

			if (row == col)
				identityMat[row * size + col] = 1.0;
			else
				identityMat[row * size + col] = 0.0;
		}
	}

	return identityMat;
}

int max(int a, int b)
{
	return a >= b ? a : b;
}

double* Conv2(double* image, int imageRow, int imageCol, double* kernel, int kernelRow, int kernelCol, char* type = "same")
{
	double* outMat;
	int outRow, outCol, edgeRows, edgeCols;

	if (!strcmp(type, "full"))
	{
		outRow = imageRow + kernelRow - 1;
		outCol = imageCol + kernelCol - 1;
		edgeRows = kernelRow - 1;
		edgeCols = kernelCol - 1;
	}
	else if (!strcmp(type, "same"))
	{
		outRow = imageRow;
		outCol = imageCol;
		edgeRows = (kernelRow - 1)/2;
		edgeCols = (kernelCol - 1)/2;
	}
	else if (!strcmp(type, "valid"))
	{
		outRow = imageRow - kernelRow + 1;
		outCol = imageCol - kernelCol + 1;
		edgeRows = edgeCols = 0;	
	}
	else
	{
		return ERROR;
	}
	
	outMat = new double[outRow * outCol];

	int iImage, iKernel, jImage, jKernel;
	double sum = 0;
	for (int i = 0; i < outRow; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			sum = 0;

			iKernel = kernelRow - 1 - max(0, edgeRows - i);
			iImage = max(0, i - edgeRows);
			for (; (iKernel>=0) && (iImage < imageRow); iKernel--, iImage++)
			{
				jKernel = kernelCol - 1 - max(0, edgeCols - j);
				jImage = max(0, j - edgeCols);

				for (; (jKernel >= 0) && (jImage < imageCol); jKernel--, jImage++)
					sum += image[imageCol * iImage + jImage] * kernel[kernelCol * iKernel + jKernel];
			}
			outMat[i * outCol + j] = sum;
		}
	}

	return outMat;
}


double* Conv2(Cell image, Cell kernel, char* type)
{
	double* outMat;
	int outRow, outCol, edgeRows, edgeCols;

	if (!strcmp(type, "full"))
	{
		outRow = image.rows + kernel.rows - 1;
		outCol = image.cols + kernel.cols - 1;
		edgeRows = kernel.rows - 1;
		edgeCols = kernel.cols - 1;
	}
	else if (!strcmp(type, "same"))
	{
		outRow = image.rows;
		outCol = image.cols;
		edgeRows = (kernel.rows - 1) / 2;
		edgeCols = (kernel.cols - 1) / 2;
	}
	else if (!strcmp(type, "valid"))
	{
		outRow = image.rows - kernel.rows + 1;
		outCol = image.cols - kernel.cols + 1;
		edgeRows = edgeCols = 0;
	}
	else
	{
		return (double*)-1;
	}

	outMat = new double[outRow * outCol];

	int iImage, iKernel, jImage, jKernel;
	double sum = 0;
	for (int i = 0; i < outRow; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			sum = 0;

			iKernel = kernel.rows - 1 - max(0, edgeRows - i);
			iImage = max(0, i - edgeRows);
			for (; (iKernel >= 0) && (iImage < image.rows); iKernel--, iImage++)
			{
				jKernel = kernel.cols - 1 - max(0, edgeCols - j);
				jImage = max(0, j - edgeCols);

				for (; (jKernel >= 0) && (jImage < image.cols); jKernel--, jImage++)
					sum += image.matx[image.cols * iImage + jImage] * kernel.matx[kernel.cols * iKernel + jKernel];
			}
			outMat[i * outCol + j] = sum;
		}
	}

	return outMat;
}

double* Fliplr(const double* arry, int height, int width) {

	double* returnBuffer = new double[width * height];
	
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
			returnBuffer[row * width + width - col - 1] = arry[row * width + col];
		}
	}

	return returnBuffer;
}

double* Flipud(const double* arry, int height, int width) {

	double* returnBuffer = new double[width * height];

	for (int col = 0; col < width; col++)
	{
		for (int row = 0; row < height; row++)
		{
			returnBuffer[(height - row - 1) * width + col] = arry[row * width + col];
		}
	}

	return returnBuffer;
}

// mat1H and mat2H are must equie
Matrix* MatrixColExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W) {

	int extWidth = mat1W + mat2W;

	Matrix* extMatrix = new Matrix;
	extMatrix->CreateMatrix(mat1H, extWidth, 1);

	for (int row = 0; row < mat1H; row++)
	{
		for (int col = 0; col < mat1W; col++)
		{
			extMatrix->mat[row * extWidth + col] = mat1[row * mat1W + col];
		}

		for (int col = mat1W - 1; col < extWidth; col++)
		{
			extMatrix->mat[row * extWidth + col] = mat2[row * mat2W + (col - (mat1W - 1))];
		}
	}

	return extMatrix;
}

// mat1W and mat2W are must equie
Matrix* MatrixRowExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W) {

	int extHeight = mat1H + mat2H;

	Matrix* extMatrix = new Matrix;
	extMatrix->CreateMatrix(extHeight, mat1W, 1);

	for (int col = 0; col < mat1W; col++)
	{
		for (int row = 0; row < mat1H; row++)
		{
			extMatrix->mat[row * mat1W + col] = mat1[row * mat1W + col];
		}
		for (int row = mat1H - 1; row < extHeight; row++)
		{
			extMatrix->mat[row * mat1W + col] = mat2[(row - (mat1H - 1)) * mat1W + col];
		}
	}

	return extMatrix;
}

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex) {

	int cutHeight = rowEndIndex - rowStartIndex + 1;
	int cutWidth = colEndIndex - colStartIndex + 1;
	
	Matrix* cutMatrix = new Matrix;
	cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);

	for (int row = rowStartIndex; row <= rowEndIndex; row++)
	{
		for (int col = colStartIndex; col <= colEndIndex; col++)
		{
			cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
		}
	}

	return cutMatrix;
}

Matrix* Upsample2df(const Matrix* h, int power) {

	int height = pow(2, power) * h->height;
	int width = pow(2, power) * h->width;

	Matrix* ho = new Matrix;
	ho->CreateMatrix(height, width, 1);

	int step = pow(2, power);
	int hStep = 0;

	for (int row = 0; row < height; row+= step)
	{
		for (int col = 0; col < width; col+= step)
		{
			ho->mat[row * ho->width + col] = h->mat[hStep];
			hStep++;
		}
	}

	return ho;
}

Matrix* Windowing(double* x, int lenghtX, int L) {

	int N = lenghtX;
	
	Matrix* y = new Matrix;
	y->height = N;
	y->width = L;
	y->depth = 1;
	y->mat = zeros(L, N);

	int T = N / L;
	double* g = zeros(2 * T, 1);

	int n = 0;
	for (int j = 0; j < 2 * T; j++)
	{
		n = -1 * T / 2 + j;
		g[j] = MeyerWind(n / T);
	}

	int index = 0;
	for (int j = 0; j < L; j++)
	{
		index = 0;
		for (int k = -1 * T / 2; k <= 1.5 * T - 1; k++)
		{
			double a = (int)(k + j * T) % N;
			int in_sig = floor(a);
			y->mat[in_sig * y->width + j] = g[index] * x[in_sig];
			index++;
		}
	}

}

double MeyerWind(double x) {

	double y = 0.0;
	
	if ((- 1 / 3 + 1 / 2 < x) && (x < 1 / 3 + 1 / 2))
		y = 1.0;

	else if (((1 / 3 + 1 / 2 <= x) && (x <= 2 / 3 + 1 / 2)) || ((-2 / 3 + 1 / 2 <= x) && (x <= 1 / 3 + 1 / 2))) {

		double w = 3 * abs(x - 1 / 2) - 1;
		double z = pow(w, 4) * (35 - 84 * w + 70 * pow(w, 2) - 20 * pow(w, 3));
		y = pow(cos(PI / 2 * (z)), 2);
	}

	else
		y = 0.0;

	return y;
}