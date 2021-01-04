#include "Container.h"

#define LINPOS(row,col,collen) (row*collen)+col

Matrix* Atrousc(Matrix* signal, Matrix* filter, double* upMatrix)
{
	double* FArray, * SArray, * outArray, * M;
	/* FArray   - Filter coefficients
	   SArray   - Signal coefficients
	   outArray - Output coefficients
	   M        - upsampling matrix 	*/
	int SColLength, SRowLength, FColLength, FRowLength, O_SColLength, O_SRowLength;
	int SFColLength, SFRowLength;
	int n1, n2, l1, l2, k1, k2, f1, f2, kk2, kk1;
	double sum;
	int M0, M3, sM0, sM3;


	SColLength = signal->width;
	SRowLength = signal->height;
	FColLength = filter->width;
	FRowLength = filter->height;

	SFColLength = FColLength - 1;
	SFRowLength = FRowLength - 1;

	FArray = filter->mat;
	SArray = signal->mat;
	M = upMatrix;
	M0 = (int)M[0];
	M3 = (int)M[3];
	sM0 = M0 - 1;
	sM3 = M3 - 1;


	O_SColLength = SColLength - M0 * FColLength + 1;
	O_SRowLength = SRowLength - M3 * FRowLength + 1;


	Matrix* outMatrix = new Matrix;
	outMatrix->CreateMatrix(O_SRowLength, O_SColLength);
	outArray = outMatrix->mat;

	/* Convolution loop */

	for (n1 = 0; n1 < O_SRowLength; n1++) {
		for (n2 = 0; n2 < O_SColLength; n2++) {
			sum = 0;
			kk1 = n1 + sM0;;
			for (k1 = 0; k1 < FRowLength; k1++) {
				kk2 = n2 + sM3;
				for (k2 = 0; k2 < FColLength; k2++) {
					f1 = SFRowLength - k1; /* flipped index */
					f2 = SFColLength - k2;
					sum += FArray[LINPOS(f1, f2, FColLength)] * SArray[LINPOS(kk1, kk2, SColLength)];
					kk2 += M3;
				}
				kk1 += M0;
			}
			outArray[LINPOS(n1, n2, O_SColLength)] = sum;
		}
	}

	return outMatrix;
}
