#pragma once
#include "Dst.h"
#include "MatlabFuncs.h"
#include <math.h>

double* Atrousrec(Dst* y, const char* lpfilt) {

	int NLevels = y->matCount - 1;

	double* h0, h1, g0, g1;
	
	// will repairing...		atrousfilters
	// [h0,h1,g0,g1] = atrousfilters(fname);

	double* x;
	double* y1;

	x = y->mats[0];
	double* I2 = Eye(2);
	double* shift = new double[2]{ 1.0, 1.0 };
	double L = 0.0;

	for (int i = NLevels - 1; i >= 1; i--) {

		y1 = y->mats[NLevels - i];			// Matlab: y{2} <=> C: y[1]
		
		shift[0] = -1 * pow(2, (i - 1)) * shift[0] + 2.0;
		shift[1] = -1 * pow(2, (i - 1)) * shift[1] + 2.0;
		
		L = pow(2, i);

		// will repairing...		symext
		// x = atrousc(symext(x, upsample2df(g0, i), shift), g0, L * I2) + atrousc(symext(y1, upsample2df(g1, i), shift), g1, L * I2);
	}

	shift[0] = 1.0;
	shift[1] = 1.0;

	// will repairing...			conv2
	// x = conv2(symext(x, g0, shift), g0, 'valid') + conv2(symext(y{ Nlevels + 1 }, g1, shift), g1, 'valid');
}