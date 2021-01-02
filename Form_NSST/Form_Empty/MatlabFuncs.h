#pragma once
#include "Cell.h"
#include "Container.h"

// Sum defined only for dim = 3 
Matrix* Sum(Matrix* mat, int dim);

double* Eye(int size);

int max(int a, int b);

// mat2 is kernel matrix.
double* Conv2(Cell image, Cell kernel, char* type = "same");

double* Conv2(double* image, int imageRow, int imageCol, double* kernel, int kernelRow, int kernelCol, char* type);

double* Fliplr(const double* arry, int height, int width);

double* Flipud(const double* arry, int height, int width);

Matrix* MatrixColExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixRowExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex);

Matrix* Upsample2df(const Matrix* mat, int power);

Matrix* Windowing(double* x, int lenghtX, int L);

double MeyerWind(double x);