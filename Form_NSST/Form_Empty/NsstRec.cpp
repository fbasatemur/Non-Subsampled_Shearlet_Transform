#include "NsstRec.h"
#include "MatlabFuncs.h"


double* NsstRec1(Cell* dst, const char* lpfilt) {

	int level = dst->matxCount - 1;
	
	Cell* y = newCell(level);
	y[0] = dst[0];

	for (int i = 1; i <= level; i++) {
		
		y[i].matx = Sum(dst[i].matx, dst[i].getSize(), dst[i].depth, 3);
	}
	
	/// will repairing...
	return y[0].matx; // real(atrousrec(&y,lpfilt));
}