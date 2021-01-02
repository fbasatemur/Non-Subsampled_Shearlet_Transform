
#include "Container.h"
#define max(a, b) (((a) > (b)) ? (a) : (b))

// Sum defined only for dim = 3 
Matrix* Sum(Matrix* mat, int dim);

double* Eye(int size);

double* ones(int width, int height);
double* zeros(int width, int height);
double* zeros(int width, int height, int depth);


Matrix* Conv2(Matrix* image, Matrix* kernel, char* type);
Matrix* Conv2(Matrix* image, Matrix kernel, char* type);
double* Conv2(double* image, int imageRow, int imageCol, double* kernel, int kernelRow, int kernelCol, char* type);

double* Fliplr(const double* arry, int height, int width);

double* Flipud(const double* arry, int height, int width);

Matrix* MatrixColExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixRowExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex);
Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex, int rowStep, int colStep);

