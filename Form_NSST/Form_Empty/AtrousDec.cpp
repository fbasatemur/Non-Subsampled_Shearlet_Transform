#pragma once
#include "AtrousDec.h"
#include <math.h>
#include "AtrousFilters.h"
#include "MatlabFuncs.h"
#include "Atrousc.h"

Cont* AtrousDec(Matrix* image, const char* lpfilt, int level)
{
	Cont* filters = AtrousFilters(lpfilt);

	double* shift = new double[2]{ 1.0, 1.0 };

	Matrix* y0 = Conv2(symext(image, filters->mats[1], shift), filters->mats[1], "valid");
	Matrix* y1 = Conv2(symext(image, filters->mats[3], shift), filters->mats[3], "valid");

	Cont* y = new Cont(level + 1);
	y->mats[level] = y0;
	image = y0;

	double* I2 = Eye(2), *I2L;
	int L;
	for (int i = 0; i < level-1; i++)
	{
		shift[0] = pow(-2, i - 1) * shift[0] + 2;
		shift[1] = pow(-2, i - 1) * shift[1] + 2;
		L = pow(2, i);
		I2L = ScalarMatMul(I2, 2 * 2, L);

		y0 = Atrousc(symext(image, Upsample2df(filters->mats[1],i), shift), filters->mats[1], I2L);
		y1 = Atrousc(symext(image, Upsample2df(filters->mats[3],i), shift), filters->mats[3], I2L);

		y->mats[level - i + 1] = y1;
		image = y0;
	}
	y->mats[0] = image;
	
	return y;
}