#pragma once
#include <Windows.h>
#include "ShearParameters.h"
#include <math.h>
#include "NsstDec.h"
#include "Cell.h"
#include "MatlabFuncs.h"

int* NsstDec1e(BYTE* image, int width, int height, struct ShearParameters shearParam, char* laplacianPyramidFilter, int* shearFilters)
{
	int level = strlen((char*)shearParam.dcomp);
	
	//Laplacian Pyramid decomposition
	//Cell* y = AtrousDec(image, laplacianPyramidFilter, level);

	Cell* dst = newCell(level + 1);
	dst[0] = y[0];

	Cell* shearF = newCell(level);
	
	for (int i = 0; i < level; i++)
	{
		//shearF[i].matx = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]) * sqrt(shearParam.dsize[i]);
		
		for (int j = 0; j < pow(2, shearParam.dcomp[i]); j++)
			dst[i + 1].matx = Conv2(y[i + 1], shearF[i]);
	}
	
}
