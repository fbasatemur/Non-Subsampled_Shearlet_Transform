#include "NsstRec.h"
#include "MatlabFuncs.h"


double* NsstRec1(Cont* dst, const char* lpfilt) {

	int level = dst->matNums - 1;
	
	Matrix** y = new Matrix*[level];

	y[0] = &dst->mats[0];


	for (int i = 1; i <= level; i++) {
		
		y[i] = Sum(&dst->mats[i], 3);
	}
	
	/// will repairing...
	return y[0]->mat; // real(atrousrec(&y,lpfilt));
}