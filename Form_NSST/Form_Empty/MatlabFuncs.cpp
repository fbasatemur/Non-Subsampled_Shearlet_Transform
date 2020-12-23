#include <string.h>		//for strcmp
#include "MatlabFuncs.h"

double* Sum(double* mat, int imageSize, int matDepth, int dim) {

	double* retMat = new double[imageSize];
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
			for (int d = 0; d < matDepth; d++)
				total += mat[d * imageSize + i];
			retMat[i] = total;
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
