#include <string.h>		//for strcmp
#include "MatlabFuncs.h"

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

double* ones(int width, int height)
{
	double* onesMat = new double[width * height];
	memset(onesMat, 1.0, height * width * sizeof(double));

	return onesMat;
}

double* zeros(int width, int height)
{
	double* zerosMat = new double[width * height];
	memset(zerosMat, 0.0, height * width * sizeof(double));

	return zerosMat;
}

double* zeros(int width, int height, int depth)
{
	double* zerosMat = new double[width * height * depth];
	memset(zerosMat, 0.0, height * width * depth * sizeof(double));

	return zerosMat;
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


Matrix* Conv2(Matrix* image, Matrix* kernel, char* type = "same")
{
	int outRow, outCol, edgeRows, edgeCols;

	if (!strcmp(type, "full"))
	{	
		outRow = image->height+ kernel->height - 1;
		outCol = image->width + kernel->width - 1;
		edgeRows = kernel->height - 1;
		edgeCols = kernel->width - 1;
	}
	else if (!strcmp(type, "same"))
	{
		outRow = image->height;
		outCol = image->width;
		edgeRows = (kernel->height - 1) / 2;
		edgeCols = (kernel->width - 1) / 2;
	}
	else if (!strcmp(type, "valid"))
	{
		outRow = image->height - kernel->height + 1;
		outCol = image->width - kernel->width + 1;
		edgeRows = edgeCols = 0;
	}
	else
	{
		return (Matrix*)-1;
	}
	
	Matrix* outMat = new Matrix;
	outMat->CreateMatrix(outRow, outCol);
	int iImage, iKernel, jImage, jKernel;
	double sum = 0;
	for (int i = 0; i < outRow; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			sum = 0;

			iKernel = kernel->height - 1 - max(0, edgeRows - i);
			iImage = max(0, i - edgeRows);
			for (; (iKernel >= 0) && (iImage < image->height); iKernel--, iImage++)
			{
				jKernel = kernel->width - 1 - max(0, edgeCols - j);
				jImage = max(0, j - edgeCols);

				for (; (jKernel >= 0) && (jImage < image->width); jKernel--, jImage++)
					sum += image->mat[image->width * iImage + jImage] * kernel->mat[kernel->width * iKernel + jKernel];
			}
			outMat->mat[i * outCol + j] = sum;
		}
	}
	return outMat;
}

Matrix* Conv2(Matrix* image, Matrix kernel, char* type = "same")
{
	int outRow, outCol, edgeRows, edgeCols;

	if (!strcmp(type, "full"))
	{
		outRow = image->height + kernel.height - 1;
		outCol = image->width + kernel.width - 1;
		edgeRows = kernel.height - 1;
		edgeCols = kernel.width - 1;
	}
	else if (!strcmp(type, "same"))
	{
		outRow = image->height;
		outCol = image->width;
		edgeRows = (kernel.height - 1) / 2;
		edgeCols = (kernel.width - 1) / 2;
	}
	else if (!strcmp(type, "valid"))
	{
		outRow = image->height - kernel.height + 1;
		outCol = image->width - kernel.width + 1;
		edgeRows = edgeCols = 0;
	}
	else
	{
		return (Matrix*)-1;
	}

	Matrix* outMat = new Matrix;
	outMat->CreateMatrix(outRow, outCol);
	int iImage, iKernel, jImage, jKernel;
	double sum = 0;
	for (int i = 0; i < outRow; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			sum = 0;

			iKernel = kernel.height - 1 - max(0, edgeRows - i);
			iImage = max(0, i - edgeRows);
			for (; (iKernel >= 0) && (iImage < image->height); iKernel--, iImage++)
			{
				jKernel = kernel.width - 1 - max(0, edgeCols - j);
				jImage = max(0, j - edgeCols);

				for (; (jKernel >= 0) && (jImage < image->width); jKernel--, jImage++)
					sum += image->mat[image->width * iImage + jImage] * kernel.mat[kernel.width * iKernel + jKernel];
			}
			outMat->mat[i * outCol + j] = sum;
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

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex, int rowStep, int colStep) {

	int cutHeight, cutWidth;
	Matrix* cutMatrix = new Matrix;

	//-,-
	if ((rowStep < 0 && rowStartIndex > rowEndIndex) && (colStep < 0 && colStartIndex > colEndIndex)) {
		cutHeight = rowStartIndex - rowEndIndex + 1;
		cutWidth = colStartIndex - colEndIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row >= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col >= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
	}
	//-,+
	else if ((rowStep < 0 && rowStartIndex > rowEndIndex) && (colStep > 0 && colStartIndex < colEndIndex)) {
		cutHeight = rowStartIndex - rowEndIndex + 1;
		cutWidth = colEndIndex - colStartIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row >= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col <= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
	}
	//+,-
	else if ((rowStep > 0 && rowStartIndex < rowEndIndex) && (colStep < 0 && colStartIndex > colEndIndex)) {
		cutHeight = rowEndIndex - rowStartIndex + 1;
		cutWidth = colStartIndex - colEndIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row <= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col >= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
	}
	//+,+
	else if ((rowStep > 0 && rowStartIndex < rowEndIndex) && (colStep > 0 && colStartIndex < colEndIndex)){
		cutHeight = rowEndIndex - rowStartIndex + 1;
		cutWidth = colEndIndex - colStartIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row <= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col <= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
	}
	return cutMatrix;
}

