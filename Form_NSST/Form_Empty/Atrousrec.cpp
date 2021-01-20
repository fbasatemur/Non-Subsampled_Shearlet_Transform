#pragma once
#include "AtrousRec.h"
#include "MatlabFuncs.h"
#include "AtrousFilters.h"
#include <math.h>
#include "Atrousc.h"

Matrix* AtrousRec(Cont* y, const char* lpfilt) {

	int NLevels = y->matNums - 1;

	// [g0, h0, g1, h1] = atrousfilters(fname);
	Cont* ret = AtrousFilters(lpfilt);
	Matrix* g0 = ret->mats[0];
	Matrix* h0 = ret->mats[1];
	Matrix* g1 = ret->mats[2];
	Matrix* h1 = ret->mats[3];


	Matrix* x;
	Matrix* y1;

	x = y->mats[0];

	Matrix* I2 = EyeMatrix(2);
	double* shift = new double[2]{ 1.0, 1.0 };
	double L = 0.0;

	for (int i = NLevels - 1; i >= 1; i--) {

		y1 = y->mats[NLevels - i];			// Matlab: y{2} <=> C: y[1]

		shift[0] = -1 * pow(2, (i - 1)) + 2.0;
		shift[1] = -1 * pow(2, (i - 1)) + 2.0;

		L = pow(2, i);

		
		// x = atrousc(symext(x, upsample2df(g0, i), shift), g0, L * I2) + atrousc(symext(y1, upsample2df(g1, i), shift), g1, L * I2);

		Matrix* mult = *I2 * L;

		Matrix* up1 = Upsample2df(g0, i);
		Matrix* symext1 = symext(x, up1, shift);
		Matrix* atrousc1 = Atrousc(symext1, g0, mult->mat);

		Matrix* up2 = Upsample2df(g1, i);
		Matrix* symext2 = symext(y1, up2, shift);
		Matrix* atrousc2 = Atrousc(symext2, g1, mult->mat);

		x = *atrousc1 + *atrousc2;

	}

	shift[0] = 1.0;
	shift[1] = 1.0;

	Matrix* symetxX = symext(x, g0, shift);
	Matrix* symetxY = symext(y->mats[NLevels], g1, shift);

	Matrix* conv2X = Conv2(symetxX, g0, "valid");
	Matrix* conv2Y = Conv2(symetxY, g1, "valid");
	x = *conv2X + *conv2Y;

	return x;
}