#pragma once
#include "Container.h"
#include "MatlabFuncs.h"
#include "AtrousFilters.h"
#include <math.h>

Matrix* Atrousrec(Cont* y, const char* lpfilt);