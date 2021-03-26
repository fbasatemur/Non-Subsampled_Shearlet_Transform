#pragma once
#include "ShearParameters.h"
#include "Container.h"

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
/// <param name="filters : ">
///		
/// </param>
/// <param name="shearFilterMyer : "></param>
/// <returns></returns>
Cont* NsstDec1e(Matrix* image, const ShearParameters& shearParam, Cont* filters, Matrix** shearFilterMyer);