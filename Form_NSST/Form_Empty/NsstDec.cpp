#pragma once
#include <math.h>
#include "NsstDec.h"
#include "MatlabFuncs.h"
#include "NSSTFuncs.h"



/// <summary>
///		This function computes the(local) nonsubsampled shearlet transform as given
/// </summary>
/// <param name="image : "> 
///		input image 
/// </param>
/// <param name="shearParam : ">
///		shear_parameters.dcomp - a vector such that.dcomp(i) indicates that the
/// 	ith decomposition level has 2 ^ decomp(i)
/// 	directions.The length of the vector plus 1 is
/// 	total the number of decompostions.
/// </param>
/// <param name="lpfilt : ">
///		lpfilt is the filter to be used for the Laplacian Pyramid / ATrous decomposition using the codes
/// </param>
/// <returns>
///		the cell array containing the discrete shearlet tranform coefficients
/// </returns>
Cont* NsstDec1e(Matrix* image, const ShearParameters& shearParam, const char* lpfilt)
{
	int level = shearParam.dcompSize;
	
	// Laplacian Pyramid decomposition	//NSLP -- cok olceklilik
	// Low and High Frequences Coefficients
	Cont* y = AtrousDec(image, lpfilt, level);

	Cont* dst = new Cont(level + 1);
	dst->mats[0] = y->mats[0];

	Cont* shearF = new Cont(level);

	for (int i = 0; i < level; i++)
	{
		int size = (int)pow(2, shearParam.dcomp[i]);	// goruntuye uygulanacak yon sayisi belirlenir.
		dst->CreateCells(i+1, size);
		shearF->CreateCells(i,size);
		
		Matrix* temp = ShearingFiltersMyer(shearParam.dsize[i], shearParam.dcomp[i]);

		for (int k = 0; k < size; k++) {
			shearF->mats[i][k] = *ScalarMatMul(temp[k], sqrt(shearParam.dsize[i]));
			dst->mats[i + 1][k] = *Conv2(y->mats[i + 1], &shearF->mats[i][k], "same");	 // Cok yonluluk 
		}	

		dst->mats[i + 1]->depth = size;

		delete[] temp;
	}

	delete shearF;

	int depth;
	for (size_t cell = 1; cell < 5; cell++) {

		depth = y->mats[cell]->depth;
		for (size_t d = 0; d < depth; d++)
			delete[] y->mats[cell][d].mat;
	}

	return dst;
}
