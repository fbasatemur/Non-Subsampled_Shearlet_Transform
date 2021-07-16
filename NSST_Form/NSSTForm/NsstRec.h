#pragma once
#include "Container.h"

/// <summary>
///		This function performs the inverse(local) nonsubsampled shearlet transform as given
/// </summary>
/// <param name="dst">
///		dst - the nonsubsampled shearlet coefficients
/// </param>
/// <param name="lpfilt">
///		lpfilt - the filter to be used for the Laplacian Pyramid / ATrous decomposition using the codes
/// </param>
/// <returns>
///		the reconstructed image
/// </returns>
Matrix* NsstRec1(Cont* dst, const Cont* filters);
