#pragma once
#include <Windows.h>
#include "ShearParameters.h"
#include <math.h>
#include "NsstDec.h"
#include "Dst.h"

int* NsstDec1e(BYTE* image, int width, int height, struct ShearParameters shearParam, char* laplacianPyramidFilter, int* shearFilters)
{
	int level = strlen((char*)shearParam.dcomp);
	
	//Laplacian Pyramid decomposition
	Dst* y = AtrousDec(image, laplacianPyramidFilter, level);

	Dst dst(level + 1, 1, level + 1);
	dst.mats[0] = y->mats[0];

	Dst shearF(level, 1, level);

	for (int i = 0; i < level; i++)
	{
		shearF.mats[i] = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]); // .* sqrt(shearParam.dsize[i]);
		for (int j = 0; j < pow(2, shearParam.dcomp[i]); j++)
		{
			//dst[i + 1] = Conv2(y[i + 1], shearF[i]);
		}
	}


	
}


Dst* AtrousDec(BYTE* image, char* laplacianPyramidFilter, int pyramidLevel)
{
}


double* ShearingFiltersMyer(int dsize, int dcomp)
{

}