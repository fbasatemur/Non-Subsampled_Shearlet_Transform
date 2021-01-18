#include "Container.h"
#define max(a, b) (((a) > (b)) ? (a) : (b))

// Sum defined only for dim = 3 
Matrix* Sum(Matrix* mat, int dim);

double* Eye(int size);
Matrix* EyeMatrix(int size);

double* ones(int width, int height);
//double* zeros(int width, int height);
double* zeros(int width, int height, int depth = 1);

Matrix* Conv2(Matrix* image, Matrix* kernel, char* type);
Matrix* Conv2(Matrix* image, Matrix kernel, char* type);
double* Conv2(double* image, int imageRow, int imageCol, double* kernel, int kernelRow, int kernelCol, char* type);

double* Fliplr(const double* arry, int height, int width);

double* Flipud(const double* arry, int height, int width);

Matrix* MatrixColExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixRowExtend(double* mat1, int mat1H, int mat1W, double* mat2, int mat2H, int mat2W);

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex);

Matrix* MatrixCut(const double* mat, int height, int width, int rowStartIndex, int rowEndIndex, int colStartIndex, int colEndIndex, int rowStep, int colStep);

Matrix* Upsample2df(const Matrix* mat, int power);

Matrix* Windowing(double* x, int lenghtX, int L);

double	MeyerWind(double x);

double* ScalarMatMul(double* mat, int matSize, double scalarValue);

double* RDivide(double* mat, double* rMat, int size);

double* MatrixMultiplication(double* m1, int row1, int col1, double* m2, int row2, int col2);

double* Linspace(int init, int finish, int N);

Matrix* GenXYCoordinates(int n);

Matrix	RecFromPol(Matrix* l, int n, Matrix* gen);

Matrix* ShearingFiltersMyer(int dcomp, int dsize);

Matrix* symext(Matrix* x, Matrix* h, double* shift);

double* AvgPol(int L, double* x1, double* y1, double* x2, double* y2);

double* FFTShift2D(double* input, int width, int height);