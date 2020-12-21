#pragma once
#include "ShearParameters.h"
#include <Windows.h>

int* NsstDec1e(BYTE* image, int width, int height, struct ShearParameters shearParam, char* laplacianPyramidFilter, int* shearFilters);