#include "Container.h"
#define max(a, b) (((a) > (b)) ? (a) : (b))


float realmod(float x, float y);

// defined for only dim = 3
Matrix* Sum(Matrix* mat, int dim);

float* Eye(int size);

Matrix* EyeMatrix(int size);

float* ones(int width, int height);

float* zeros(int width, int height, int depth = 1);

float* Conv2(float* image, int imageRow, int imageCol, float* kernel, int kernelRow, int kernelCol, char* type = "same");

Matrix* Conv2(Matrix* image, Matrix* kernel, char* type = "same");

float* Fliplr(const float* arry, int height, int width);

float* Flipud(const float* arry, int height, int width);

Matrix* MatrixColExtend(float* mat1, int mat1H, int mat1W, float* mat2, int mat2H, int mat2W);

Matrix* MatrixRowExtend(float* mat1, int mat1H, int mat1W, float* mat2, int mat2H, int mat2W);

Matrix* MatrixCut(const float* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex);

Matrix* MatrixCut(const float* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex, int rowStep, int colStep);

float* ScalarMatMul(float* mat, int matSize, float scalarValue);

Matrix*	ScalarMatMul(Matrix& matx, float scalarValue);

float* RDivide(float* mat, float* rMat, int size);

float* MatrixMultiplication(float* m1, int row1, int col1, float* m2, int row2, int col2);

float* FFTShift2D(float* input, int width, int height);