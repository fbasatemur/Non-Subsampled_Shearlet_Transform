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

	Matrix* temp;
	for (int i = 0; i < level; i++)
	{
		
		temp = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]);
		//BURAYA KADAR SORUN YOK.
		temp->mat = ScalarMatMul(temp->mat, temp->GetSize2D(), sqrt(shearParam.dsize[i]));
		shearF->mats[i] = temp;

		for (int k = 0; k < pow(2, shearParam.dcomp[i]); k++)
			dst->mats[i + 1] = Conv2(y->mats[i+1], shearF->mats[i][k], "same");
	}
	
	return dst;
}
