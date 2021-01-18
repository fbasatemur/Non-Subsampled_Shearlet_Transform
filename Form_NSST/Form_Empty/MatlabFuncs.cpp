#include <string.h>		//for strcmp
#include "MatlabFuncs.h"
#include <cmath>
#include "Process.h"

# define PI           3.14159265358979323846
#define ERROR (double*)-1

inline double realmod(double x, double y)
{
	double result = fmod(x, y);
	return result >= 0 ? result : result + y;
}

Matrix* Sum(Matrix* mat, int dim) {

	int imageSize = mat->GetSize2D();

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

Matrix* EyeMatrix(int size) {

	Matrix* identityMat = new Matrix(size, size, 1);

	for (int row = 0; row < size; row++) {
		for (int col = 0; col < size; col++) {

			if (row == col)
				identityMat->mat[row * size + col] = 1.0;
			else
				identityMat->mat[row * size + col] = 0.0;
		}
	}

	return identityMat;
}

double* ones(int width, int height)
{
	double* onesMat = new double[width * height];
	for (int i = 0; i < width * height; i++)
		onesMat[i] = 1.0;

	return onesMat;
}

//double* zeros(int width, int height)
//{
//	double* zerosMat = new double[width * height];
//	for (int i = 0; i < width * height; i++)
//		zerosMat[i] = 0.0;
//
//	return zerosMat;
//}

double* zeros(int width, int height, int depth)
{
	//double* zerosMat = new double[width * height * depth];
	//for (int i = 0; i < height * width * depth; i++)
	//	zerosMat[i] = 0.0;

	return (double*)calloc(width * height * depth, sizeof(double));
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

		for (int col = mat1W; col < extWidth; col++)
		{
			extMatrix->mat[row * extWidth + col] = mat2[row * mat2W + (col - mat1W)];
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
		for (int row = mat1H; row < extHeight; row++)
		{
			extMatrix->mat[row * mat1W + col] = mat2[(row - (mat1H)) * mat1W + col];
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
			cutMatrix->mat[(row - rowStartIndex) * cutWidth + (col - colStartIndex)] = mat[row * width + col];
		}
	}

	return cutMatrix;
}


Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex, int rowStep, int colStep) {

	int cutHeight, cutWidth;
	Matrix* cutMatrix = new Matrix;

	//-,-
	if ((rowStep < 0 && rowStartIndex >= rowEndIndex) && (colStep < 0 && colStartIndex >= colEndIndex)) {
		cutHeight = rowStartIndex - rowEndIndex + 1;
		cutWidth = colStartIndex - colEndIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row >= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col >= colEndIndex; col += colStep)
				cutMatrix->mat[(rowStartIndex - row) * cutWidth + (colStartIndex - col)] = mat[row * width + col];
	}
	//-,+
	else if ((rowStep < 0 && rowStartIndex >= rowEndIndex) && (colStep > 0 && colStartIndex <= colEndIndex)) {
		cutHeight = rowStartIndex - rowEndIndex + 1;
		cutWidth = colEndIndex - colStartIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row >= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col <= colEndIndex; col += colStep)
				cutMatrix->mat[(rowStartIndex - row) * cutWidth + col] = mat[row * width + col];
	}
	//+,-
	else if ((rowStep > 0 && rowStartIndex <= rowEndIndex) && (colStep < 0 && colStartIndex >= colEndIndex)) {
		cutHeight = rowEndIndex - rowStartIndex + 1;
		cutWidth = colStartIndex - colEndIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row <= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col >= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + (colStartIndex - col)] = mat[row * width + col];
	}
	//+,+
	else if ((rowStep > 0 && rowStartIndex <= rowEndIndex) && (colStep > 0 && colStartIndex <= colEndIndex)){
		cutHeight = rowEndIndex - rowStartIndex + 1;
		cutWidth = colEndIndex - colStartIndex + 1;
		cutMatrix->CreateMatrix(cutHeight, cutWidth, 1);
		for (int row = rowStartIndex; row <= rowEndIndex; row += rowStep)
			for (int col = colStartIndex; col <= colEndIndex; col += colStep)
				cutMatrix->mat[row * cutWidth + col] = mat[row * width + col];
	}
	return cutMatrix;
}


Matrix* Upsample2df(const Matrix* h, int power) {

	int height = pow(2, power) * h->height;
	int width = pow(2, power) * h->width;

	Matrix* ho = new Matrix;
	ho->mat = zeros(width, height);
	ho->width = width; ho->height = height; ho->depth = 1;

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

	double T = N / L;
	double* g = zeros(2 * T, 1);

	double n = 0;
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
			int in_sig = floor(realmod((int)(k + j * T), N));
			y->mat[in_sig * y->width + j] = g[index] * x[in_sig];
			index++;
		}
	}

	return y;
}

double MeyerWind(double x) {

	double y = 0.0;
	
	if ((- 1.0 / 3.0 + 1.0 / 2.0 < x) && (x < 1.0 / 3.0 + 1.0 / 2.0))
		y = 1.0;

	else if (((1.0 / 3.0 + 1.0 / 2.0 <= x) && (x <= 2.0 / 3.0 + 1.0 / 2.0)) || ((-2.0 / 3.0 + 1.0 / 2.0 <= x) && (x <= 1.0 / 3.0 + 1.0 / 2.0))) {

		double w = 3.0 * abs(x - 1.0 / 2.0) - 1.0;
		double z = pow(w, 4) * (35 - 84 * w + 70 * pow(w, 2) - 20 * pow(w, 3));
		y = pow(cos(PI / 2.0 * (z)), 2);
	}

	else
		y = 0.0;

	return y;
}

double* ScalarMatMul(double* mat, int matSize, double scalarValue)
{
	double* mat2 = new double[matSize]();
	for (int i = 0; i < matSize; i++)
		mat2[i] = mat[i] *scalarValue;
	return mat2;
}

double* RDivide(double* mat, double* rMat, int size)
{
	double* retMat = new double[size];
	for (int i = 0; i < size; i++)
		retMat[i] = mat[i] / rMat[i];

	return retMat;
}


double* MatrixMultiplication(double* m1, int row1, int col1, double* m2, int row2, int col2)
{
	if (col1 != row2)
		exit(1);

	double sum = 0;
	double* m3 = new double[row1 * col2];
	for (int i = 0; i < row1; i++) {
		for (int j = 0; j < col2; j++) {
			sum = 0;
			for (int k = 0; k < col1; k++) {
				sum += m1[i * col1 + k] * m2[k * col2 + j];
			}
			m3[i * col2 + j] = sum;
		}
	}
	return m3;
}

double* Linspace(int d1, int d2, int N)
{
	double a1 = d1, a2 = d2;

	int n1 = N - 1;
	double* y = new double[N]();

	for (int i = 0; i <= n1; i++)
		y[i] = a1 + i * (a2 - a1) / n1;

	y[0] = a1; y[n1] = a2;

	return y;
}


Matrix* GenXYCoordinates(int n)
{
	++n;
	double* x1 = zeros(n, n);
 	double* y1 = zeros(n, n);
	double* x2 = zeros(n, n);
	double* y2 = zeros(n, n);

	double* xt = zeros(1, n);
	double* m  = zeros(1, n);

	double y0 = 1.0, x0, xN, yN;
	int flag;
	for (int i = 0; i < n; i++)
	{
		x0 = i+1; 
		xN = n - x0 + 1; 
		yN = n;

		if (xN == x0)   flag = 1;
		else {
			m[i] = (yN - y0) / (xN - x0);	
			flag = 0;
		}

		xt = Linspace(x0, xN, n);

		int index;
		for (int j = 0; j < n; j++)
		{
			index = i * n + j;
			if (flag == 0)
			{
				y1[index] = round(m[i] * (xt[j] - x0) + y0);
				x1[index] = round(xt[j]);
				x2[index] = y1[index];
				y2[index] = x1[index];
			}
			else
			{
				x1[index] = (n - 1) / 2 + 1;
				y2[index] = (n - 1) / 2 + 1;
				y1[index] = j + 1;
				x2[index] = j + 1;
			}
		}

	}

	--n;
	double* x1n = zeros(n, n);
	double* y1n = zeros(n, n);
	double* x2n = zeros(n, n);
	double* y2n = zeros(n, n);

	Matrix* x1Cut = MatrixCut(x1, n + 1, n + 1, 0, n-1, 0, n - 1);
	Matrix* y1Cut = MatrixCut(y1, n + 1, n + 1, 0, n-1, 0, n - 1);
	Matrix* x2Cut = MatrixCut(x2, n + 1, n + 1, 1, n, 0, n - 1);
	Matrix* y2Cut = MatrixCut(y2, n + 1, n + 1, 1, n, 0, n - 1);

	int index;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			index = i * n + j;
			x1n[index] = x1Cut->mat[index];
			y1n[index] = y1Cut->mat[index];
			x2n[index] = x2Cut->mat[index];
			y2n[index] = y2Cut->mat[index];
		}

	//correct for portion outside boundry
	x1n = Flipud(x1n, n, n);
	y2n[(n-1)*n] = n;

	//return [x1n,y1n,x2n,y2n,D]
	Matrix* ret = new Matrix[5];
	ret[0].mat = x1n; ret[0].width = n; ret[0].height = n; ret[0].depth = 1; 
	ret[1].mat = y1n; ret[1].width = n; ret[1].height = n; ret[0].depth = 1; 
	ret[2].mat = x2n; ret[2].width = n; ret[2].height = n; ret[0].depth = 1; 
	ret[3].mat = y2n; ret[3].width = n; ret[3].height = n; ret[0].depth = 1; 

	ret[4].mat = AvgPol(n, x1n, y1n, x2n, y2n);
	ret[4].width = n; ret[4].height = n;

	return ret;
}

double* AvgPol(int L, double* x1, double* y1, double* x2, double* y2)
{
	double* D = zeros(L, L);

	//int offset;
	for (int i = 0; i < L; i++)
		for (int j = 0; j < L; j++)
			D[(int)((y1[i * L + j]-1) * L + (x1[i * L + j]-1))]++;
			
	for (int i = 0; i < L; i++)
		for (int j = 0; j < L; j++)
			D[(int)((y2[i * L + j]-1) * L + (x2[i * L + j]-1))]++;
		
	return D;
}

Matrix RecFromPol(Matrix* l, int n, Matrix* gen)
{
	int offset;
	double* C = zeros(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			offset = (gen[1].mat[i * n + j]-1) * n + (gen[0].mat[i * n + j]-1);
			C[offset] += l->mat[i * n + j];

			offset = (gen[3].mat[i * n + j]-1) * n + (gen[2].mat[i * n + j]-1);
			C[offset] += l->mat[(i + n) * n + j];
		}

	C = RDivide(C, gen[4].mat, n * n);

	Matrix matC;
	matC.mat = C;   matC.width = n;    matC.height = n;		matC.depth = 1;

	return matC;
}


Matrix* ShearingFiltersMyer(int n, int level)
{
	//gen : [x11, y11, x12, y12, F1]
	Matrix* gen = GenXYCoordinates(n);

	Matrix* wf = Windowing(ones(2 * n, 1), 2 * n, pow(2, level));

	int size = pow(2, level);


	Matrix* wS = new Matrix[size];
	for (int i = 0; i < size; i++)
	{
		wS[i].width = n;
		wS[i].height = n;
		wS[i].depth = 1;
		wS[i].mat = zeros(n, n);
	}

	double* one = ones(n, 1);
	Matrix* temp;
	double* fftS1;
	double* inImag  = zeros(n, n);
	double* outReal = zeros(n, n);
	double* outImag = zeros(n, n);
	
	for (int i = 0; i < size; i++)
	{
		temp = MatrixCut(wf->mat, wf->height, wf->width, 0, wf->height-1, i, i);
		temp->mat = MatrixMultiplication(temp->mat, temp->height, temp->width, one, 1, n);

		wS[i] = RecFromPol(temp, n, gen);

		// w_s(:, : , k) = real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);
		
		fftS1 = FFTShift2D(wS[i].mat, wS[i].width, wS[i].height);				// fftshift(w_s(:, : , k))

		IFFT2D(outReal, outImag, fftS1, inImag, wS[i].width, wS[i].height);		// ifft2(fftshift(w_s(:, : , k)))

		wS[i].mat = FFTShift2D(outReal, wS[i].width, wS[i].height);				// real(fftshift(ifft2(fftshift(w_s(:, : , k)))))
		
		wS[i] /= (double)sqrt(n);												// real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);

	}

	return wS;
}


Matrix* symext(Matrix* x, Matrix* h, double* shift)
{
	int m = x->width;
	int n = x->height;
	int p = h->width;
	int q = h->height;

	double p2 = floor(p / 2);
	double q2 = floor(q / 2);

	double s1 = shift[0];
	double s2 = shift[1];

	double ss = p2 - s1 + 1;
	double rr = q2 - s2 + 1;


	Matrix* yT, *temp, *extentedMatrix;

	//[fliplr(x(:,1:ss)) x  x(:,n  :-1: n-p-s1+1)]
	temp			= MatrixCut(x->mat, x->height, x->width, 0, x->height-1, 0, ss-1); // x(:,1:ss)
	temp->mat		= Fliplr(temp->mat, temp->height, temp->width);// fliplr(x(:,1:ss))
	extentedMatrix	= MatrixColExtend(temp->mat, temp->height, temp->width, x->mat, x->height, x->width);
	temp			= MatrixCut(x->mat, x->height, x->width, 0, x->height-1, n-1, n - p - s1 + 1 - 1, 1, -1); // x(:, n : -1 : n - p - s1 + 1)
	yT				= MatrixColExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);

	//[flipud(yT(1:rr, : )); yT;  yT(m  :-1 : m - q - s2 + 1, : )]
	temp			= MatrixCut(yT->mat, yT->height, yT->width, 0, rr-1, 0, yT->width-1); //yT(1:rr, : )
	temp->mat		= Flipud(temp->mat, temp->height, temp->width);		//flipud(yT(1:rr, : ))
	extentedMatrix	= MatrixRowExtend(temp->mat, temp->height, temp->width, yT->mat, yT->height, yT->width);
	temp			= MatrixCut(yT->mat, yT->height, yT->width, m-1, m - q - s2 + 1 -1, 0, yT->width-1, -1, 1);	//yT(m  :-1 : m - q - s2 + 1, : )
	yT				= MatrixRowExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);

	// yT(1:m+p-1 ,1:n+q-1)
	yT = MatrixCut(yT->mat, yT->height, yT->width, 0, m + p - 1-1, 0, n + q - 1-1);

	return yT;
}

void circshift(double* output, double* input, int rows, int cols, int yshift, int xshift)
{
	for (int r = 0; r < rows; r++) {

		int newR = (r + yshift) % rows;
		
		if (newR < 0) 
			newR = rows + newR;
		
		for (int c = 0; c < cols; c++) {

			int newC = (c + xshift) % cols;

			if (newC < 0) 
				newC = cols + newC;

			output[newR * cols + newC] = input[r * cols + c];
		}
	}
}

double* FFTShift2D(double* input, int width, int height){

	double* temp = new double[width * height];

	// x = circshift(x,floor(size(x)/2));
	circshift(temp, input, height, width, floor(height / 2), floor(width / 2));

	return temp;
}