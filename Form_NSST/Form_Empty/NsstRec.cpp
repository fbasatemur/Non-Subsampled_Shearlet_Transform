#include "NsstRec.h"

double* NsstRec1(double** dst, const char* lpfilt, int dstLenght) {

	int level = dstLenght - 1;
	double** y = new double* [dstLenght];
	y[0] = dst[0];

	/*for (int i = 1; i <= level; i++) {
		y[i] = dst[i];
	}

	return x;*/
}