#include "symext.h"
#include <cmath> //for floor
#include "MatlabFuncs.h"

Matrix* symetx(Matrix* x, Matrix* h, double* shift)
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


	Matrix* yT = new Matrix;
	Matrix* temp,* extentedMatrix;

	//[fliplr(x(:,1:ss)) x  x(:,n  :-1: n-p-s1+1)]
	temp			= MatrixCut(x->mat, x->height, x->width, 0, x->height, 0, ss); // x(:,1:ss)
	temp->mat		= Fliplr(temp->mat, x->height, x->width);// fliplr(x(:,1:ss))
	extentedMatrix  = MatrixColExtend(temp->mat, temp->height, temp->width, x->mat, x->height, x->width);
	temp			= MatrixCut(x->mat, x->height, x->width, 0, x->height, n, n - p - s1 + 1, 1, -1); // x(:, n : -1 : n - p - s1 + 1)
	yT				= MatrixColExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);

	//[flipud(yT(1:rr, : )); yT;  yT(m  :-1 : m - q - s2 + 1, : )]
	temp			= MatrixCut(yT->mat, yT->height, yT->width, 1, rr, 0, yT->height); //yT(1:rr, : )
	temp->mat		= Flipud(temp->mat, temp->height, temp->width);		//flipud(yT(1:rr, : ))
	extentedMatrix	= MatrixRowExtend(temp->mat, temp->height, temp->width, yT->mat, yT->height, yT->width);
	temp			= MatrixCut(yT->mat, yT->height, yT->width, m, m - q - s2 + 1, 0, yT->height, -1, 1);	//yT(m  :-1 : m - q - s2 + 1, : )
	yT				= MatrixRowExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);

	// yT(1:m+p-1 ,1:n+q-1)
	yT = MatrixCut(yT->mat, yT->height, yT->width, 0, m + p - 1, 0, n + q - 1);

	return yT;
}
