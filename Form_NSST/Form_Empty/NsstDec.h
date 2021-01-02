#pragma once
#include "ShearParameters.h"
#include <Windows.h>
#include "Container.h"

Cont* NsstDec1e(Matrix* image, struct ShearParameters shearParam, const char* lpfilt);