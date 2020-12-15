#include "NsstRec.h"
#include "MatlabFuncs.h"


double* NsstRec1(Dst* dst, const char* lpfilt) {

	int level = dst->matCount - 1;
	
	Dst y(dst->cols, dst->rows, dst->matCount);
	
	for (int i = 0; i < dst->matCount; i++)
		y.setMatDepth(i, 1);

	y.mats[0] = dst->mats[0];

	for (int i = 1; i <= level; i++) {
		
		y.mats[i] = Sum(dst->mats[i], dst->getSize(), dst->matsDepth[i], 3);
	}
	
	/// will repairing...
	return y.mats[0]; // real(atrousrec(&y,lpfilt));
}