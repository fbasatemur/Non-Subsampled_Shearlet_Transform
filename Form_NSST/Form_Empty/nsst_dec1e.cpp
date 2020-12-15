#include <Windows.h>
#include "ShearParameters.h"
#include <math.h>

int* NsstDec1e(BYTE* image, int width, int height, struct ShearParameters shearParam, char* laplacianPyramidFilter, int* shearFilters)
{
	int pyramidLevel = strlen((char*)shearParam.dcomp);
	
	//Laplacian Pyramid  decomposition
	double* y = AtrousDec(image, laplacianPyramidFilter, pyramidLevel);

	double** dst = (double**)malloc((pyramidLevel+1) * sizeof(double*));
	dst[0] = y;

	double** shearF = (double**)malloc(pyramidLevel * sizeof(double*));

	for (int i = 0; i < pyramidLevel; i++)
	{
		shearF[i] = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]); // .* sqrt(shearParam.dsize[i]);
		for (int j = 0; j < pow(2, shearParam.dcomp[i]); j++)
		{
			//dst[i + 1] = Conv2(y[i + 1], shearF[i]);
		}
	}


	
}


double* AtrousDec(BYTE* image, char* laplacianPyramidFilter, int pyramidLevel)
{
}


double* ShearingFiltersMyer(int dsize, int dcomp)
{

}