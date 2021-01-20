#pragma once
#include <math.h>
#include <Windows.h>
#include "NsstDec.h"
#include "MatlabFuncs.h"
#include "AtrousDec.h"

Cont* NsstDec1e(Matrix* image, const ShearParameters& shearParam, const char* lpfilt)
{
	int level = shearParam.dcompSize;
	
	//Laplacian Pyramid decomposition
	Cont* y = AtrousDec(image, lpfilt, level);

	Cont* dst = new Cont(level + 1);
	dst->mats[0] = y->mats[0];

	Cont* shearF = new Cont(level);
	// shearF->CreateCells();

	Matrix* temp;
	for (int i = 0; i < level; i++)
	{
		int size = pow(2, shearParam.dcomp[i]);
		dst->CreateCells(i+1, size);
		shearF->CreateCells(i,size);

		temp = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]);

		for (int k = 0; k < size; k++) {
			shearF->mats[i][k] = ScalarMatMul(temp[k], sqrt(shearParam.dsize[i]));
			dst->mats[i + 1][k] = Conv2(y->mats[i + 1], shearF->mats[i][k], "same");
		}	

		dst->mats[i + 1]->depth = size;
	}
	return dst;
}
