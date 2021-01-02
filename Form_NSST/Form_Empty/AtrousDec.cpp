#pragma once
#include "AtrousDec.h"
#include <math.h>
#include "symext.h"
#include "AtrousFilters.h"
#include "MatlabFuncs.h"

Cont* AtrousDec(Matrix* image, const char* lpfilt, int level)
{
	Cont* filters = AtrousFilters(lpfilt);

	Cont* y = new Cont(level + 1);

	double* shift = new double[2]{ 1.0, 1.0 };

	Matrix* y0 = Conv2(symetx(image, filters->mats[1], shift), filters->mats[1], "valid");
	Matrix* y1 = Conv2(symetx(image, filters->mats[3], shift), filters->mats[3], "valid");

	y->mats[level] = y0;
	image = y0;

	double* I2 = Eye(2);
	int L;
	for (int i = 0; i < level-1; i++)
	{
		shift[0] = pow(-2, i - 1) * shift[0] + 2;
		shift[1] = pow(-2, i - 1) * shift[1] + 2;
		L = pow(2, i);

		//y0 = atrousc(symext(x,upsample2df(h0,i),shift),h0,I2 * L);
		//y1 = atrousc(symext(x,upsample2df(h1,i),shift),h1,I2 * L);

		y->mats[level - i + 1] = y1;
		image = y0;
	}
	y->mats[0] = image;
	
	return y;
}
