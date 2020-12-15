#include "MatlabFuncs.h"
#include "Dst.h"

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