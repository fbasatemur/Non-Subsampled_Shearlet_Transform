#include "MatlabFuncs.h"
#include "Process.h"
#include <cmath>
#define PI           3.14159265358979323846
#define LINPOS(row,col,collen) (row*collen)+col

float* Linspace(int d1, int d2, int N)
{
	float a1 = d1, a2 = d2;

	int n1 = N - 1;
	float* y = new float[N]();

	for (int i = 0; i <= n1; i++)
		y[i] = a1 + i * (a2 - a1) / n1;

	y[0] = a1; y[n1] = a2;

	return y;
}

/// <summary>
///		This function generates the matrix that contains the number of times the polar grid points go through the rectangular grid 
/// </summary>
/// <param name="L : ">
///		L is the order of the block matrix
/// </param>
/// <returns>
///		D is the common grid point values
/// </returns>
float* AvgPol(int L, float* x1, float* y1, float* x2, float* y2)
{
	float* D = zeros(L, L);

	//int offset;
	for (int i = 0; i < L; i++)
		for (int j = 0; j < L; j++)
			D[(int)((y1[i * L + j] - 1) * L + (x1[i * L + j] - 1))]++;

	for (int i = 0; i < L; i++)
		for (int j = 0; j < L; j++)
			D[(int)((y2[i * L + j] - 1) * L + (x2[i * L + j] - 1))]++;

	return D;
}

Matrix* Upsample2df(const Matrix* h, int power) {

	int height = pow(2, power) * h->height;
	int width = pow(2, power) * h->width;

	Matrix* ho = new Matrix(height, width);
	ho->mat = zeros(width, height);

	int step = pow(2, power);
	int hStep = 0;

	for (int row = 0; row < height; row += step)
	{
		for (int col = 0; col < width; col += step)
		{
			ho->mat[row * ho->width + col] = h->mat[hStep];
			hStep++;
		}
	}

	return ho;
}


/// <summary>
///		This function generates the xand y vectors that contain the i, j coordinates to extract radial slices
/// </summary>
/// <param name="n : ">
///		n is the order of the block to be used
/// </param>
/// <returns>
///		** return [x1n, y1n, x2n, y2n, D] ** <para></para>
///		x1, y1 are the i, j values that correspond to the radial slices from the endpoints 1, 1 to n, 1 through the origin <para></para>
///		x2, y2 are the i, j values that correspond to the radial slices from the endpoints 1, 1 to 1, n through the origin <para></para>
///		D is the matrix that contains the number of times the polar grid points go through the rectangular grid
/// </returns>
Matrix* GenXYCoordinates(int n)
{
	++n;
	float* x1 = zeros(n, n);
	float* y1 = zeros(n, n);
	float* x2 = zeros(n, n);
	float* y2 = zeros(n, n);

	float* xt = zeros(1, n);
	float* m = zeros(1, n);

	float y0 = 1.0, x0, xN, yN;
	int flag;
	for (int i = 0; i < n; i++)
	{
		x0 = i + 1;
		xN = n - x0 + 1;
		yN = n;

		if (xN == x0)   flag = 1;
		else {
			m[i] = (yN - y0) / (xN - x0);
			flag = 0;
		}

		xt = Linspace(x0, xN, n);

		int index;
		for (int j = 0; j < n; j++)
		{
			index = i * n + j;
			if (flag == 0)
			{
				y1[index] = round(m[i] * (xt[j] - x0) + y0);
				x1[index] = round(xt[j]);
				x2[index] = y1[index];
				y2[index] = x1[index];
			}
			else
			{
				x1[index] = (n - 1) / 2 + 1;
				y2[index] = (n - 1) / 2 + 1;
				y1[index] = j + 1;
				x2[index] = j + 1;
			}
		}

	}

	--n;
	float* x1n = zeros(n, n);
	float* y1n = zeros(n, n);
	float* x2n = zeros(n, n);
	float* y2n = zeros(n, n);

	Matrix* x1Cut = MatrixCut(x1, n + 1, n + 1, 0, n - 1, 0, n - 1);
	Matrix* y1Cut = MatrixCut(y1, n + 1, n + 1, 0, n - 1, 0, n - 1);
	Matrix* x2Cut = MatrixCut(x2, n + 1, n + 1, 1, n, 0, n - 1);
	Matrix* y2Cut = MatrixCut(y2, n + 1, n + 1, 1, n, 0, n - 1);

	int index;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			index = i * n + j;
			x1n[index] = x1Cut->mat[index];
			y1n[index] = y1Cut->mat[index];
			x2n[index] = x2Cut->mat[index];
			y2n[index] = y2Cut->mat[index];
		}

	//correct for portion outside boundry
	float* x1n_ = Flipud(x1n, n, n); delete[] x1n;
	y2n[(n - 1) * n] = n;

	// return [x1n,y1n,x2n,y2n,D]
	Matrix* ret = new Matrix[5];
	ret[0].mat = x1n_; ret[0].width = n; ret[0].height = n; ret[0].depth = 1;
	ret[1].mat = y1n; ret[1].width = n; ret[1].height = n; ret[0].depth = 1;
	ret[2].mat = x2n; ret[2].width = n; ret[2].height = n; ret[0].depth = 1;
	ret[3].mat = y2n; ret[3].width = n; ret[3].height = n; ret[0].depth = 1;

	ret[4].mat = AvgPol(n, x1n_, y1n, x2n, y2n);
	ret[4].width = n; ret[4].height = n;

	delete[] x1; delete[] y1; delete[] x2; delete[] y2; delete[] xt; delete[] m;

	return ret;
}


/// <summary>
///		This funcion re - assembles the radial slice into block.
/// </summary>
/// <param name="l : ">
///		l is the radial slice matrix
/// </param>
/// <param name="n : ">
///		n is the order of the matrix that is to be re - assembled
///		x1, y1, x2, y2 are the polar coordinates generated from function
/// </param>
/// <param name="gen : ">
///		D is the matrix containing common polar grid points
/// </param>
/// <returns>
///		C is the re - assembled block matrix
/// </returns>
Matrix* RecFromPol(Matrix* l, int n, Matrix* gen)
{
	int offset;
	float* C = zeros(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			offset = (gen[1].mat[i * n + j] - 1) * n + (gen[0].mat[i * n + j] - 1);
			C[offset] += l->mat[i * n + j];

			offset = (gen[3].mat[i * n + j] - 1) * n + (gen[2].mat[i * n + j] - 1);
			C[offset] += l->mat[(i + n) * n + j];
		}

	float* retC = RDivide(C, gen[4].mat, n * n);
	delete[] C;

	Matrix* matC = new Matrix(n, n, 1);
	matC->mat = retC;

	return matC;
}


/// <summary>
///		This function computes the Meyer window function of signal x
/// </summary>
/// <param name="x : ">
///		x - signal to be windowed
/// </param>
/// <returns>
///		y - windowed signal
/// </returns>
float MeyerWind(float x) {

	float y = 0.0;

	if ((-1.0 / 3.0 + 1.0 / 2.0 < x) && (x < 1.0 / 3.0 + 1.0 / 2.0))
		y = 1.0;

	else if (((1.0 / 3.0 + 1.0 / 2.0 <= x) && (x <= 2.0 / 3.0 + 1.0 / 2.0)) || ((-2.0 / 3.0 + 1.0 / 2.0 <= x) && (x <= 1.0 / 3.0 + 1.0 / 2.0))) {

		float w = 3.0 * abs(x - 1.0 / 2.0) - 1.0;
		float z = pow(w, 4) * (35 - 84 * w + 70 * pow(w, 2) - 20 * pow(w, 3));
		y = pow(cos(PI / 2.0 * (z)), 2);
	}

	else y = 0.0;

	return y;
}


/// <summary>
///		This function computes L band sequence of a Meyer windowed signal.
/// </summary>
/// <param name="x : ">
///		x - signal to be decomposed
/// </param>
/// <param name="lenghtX : "></param>
/// <param name="L : ">
///		L - number of bandpass filters
/// </param>
/// <returns>
///		y - windowed signal of x where second index indicates which bandpass index
/// </returns>
Matrix* Windowing(float* x, int lenghtX, int L) {

	int N = lenghtX;

	Matrix* y = new Matrix;
	y->height = N;
	y->width = L;
	y->depth = 1;
	y->mat = zeros(L, N);

	float T = N / L;
	float* g = zeros(2 * T, 1);

	float n = 0;
	for (int j = 0; j < 2 * T; j++)
	{
		n = -1 * T / 2 + j;
		g[j] = MeyerWind(n / T);
	}

	int index = 0;
	for (int j = 0; j < L; j++)
	{
		index = 0;
		for (int k = -1 * T / 2; k <= 1.5 * T - 1; k++)
		{
			int in_sig = floor(realmod((int)(k + j * T), N));
			y->mat[in_sig * y->width + j] = g[index] * x[in_sig];
			index++;
		}
	}

	delete[] g;

	return y;
}


/// <summary>
///		This function computes the directional / shearing filters using the Meyer window function.
/// </summary>
/// <param name="n : ">
///		n indicates the supports size of the directional filter is n1xn1 level indicates that the number of directions is 2 ^ level
/// </param>
/// <param name="level : ">
///		level indicates
/// </param>
/// <returns>
///		a sequence of 2D directional/shearing filters w_s where the third index determines which directional filter is used
/// </returns>
Matrix* ShearingFiltersMyer(int n, int level)
{
	//Pseudo-polar koordinat sistemi için indexlerin uretilmesi
	//gen : [x11, y11, x12, y12, F1]
	Matrix* gen = GenXYCoordinates(n);
	float* ones_ = ones(2 * n, 1);
	Matrix* wf = Windowing(ones_, 2 * n, pow(2, level));
	delete[] ones_;

	int size = pow(2, level);

	Matrix* wS = new Matrix[size];

	float* one = ones(n, 1);
	float* fftS1;
	float* inImag = zeros(n, n);
	float* outReal = zeros(n, n);
	float* outImag = zeros(n, n);

	Matrix* temp, * temp2;
	for (int i = 0; i < size; i++)
	{
		temp = MatrixCut(wf->mat, wf->height, wf->width, 0, wf->height - 1, i, i);

		temp2 = new Matrix(temp->height, temp->width);
		temp2->mat = MatrixMultiplication(temp->mat, temp->height, temp->width, one, 1, n);

		wS[i] = *RecFromPol(temp2, n, gen);

		// w_s(:, : , k) = real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);

		fftS1 = FFTShift2D(wS[i].mat, wS[i].width, wS[i].height);				// fftshift(w_s(:, : , k))

		IFFT2D(outReal, outImag, fftS1, inImag, wS[i].width, wS[i].height);		// ifft2(fftshift(w_s(:, : , k)))

		//window dizisini kartezyen koordinat sistemine donusturur.
		wS[i].mat = FFTShift2D(outReal, wS[i].width, wS[i].height);				// real(fftshift(ifft2(fftshift(w_s(:, : , k)))))

		wS[i] /= (float)sqrt(n);												// real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);

		delete temp; delete[] fftS1; delete temp2;
	}

	delete[] one; delete[] outReal; delete[] outImag; delete[] inImag;
	delete[] gen; delete wf;

	return wS;
}


/// <summary>
///		Performs symmetric extension for image x, filter h. <para></para>
///		The filter h is assumed have odd dimensions. <para></para> 
///		If the filter has horizontaland vertical symmetry, then the nonsymmetric part of conv2(h, x) has the same size of x.
/// </summary>
/// <param name="x : ">
///		mxn image
/// </param>
/// <param name="h : ">
///		2 - D filter coefficients
/// </param>
/// <param name="shift : ">
///		optional shift
/// </param>
/// <returns>
///		yT image symetrically extended(H / V symmetry)
/// </returns>
Matrix* symext(Matrix* x, Matrix* h, float* shift)
{
	int m = x->height;
	int n = x->width;
	int p = h->height;
	int q = h->width;

	float p2 = floor(p / 2);
	float q2 = floor(q / 2);

	float s1 = shift[0];
	float s2 = shift[1];

	float ss = p2 - s1 + 1;
	float rr = q2 - s2 + 1;


	// returnleri incele, parametleri sil
	Matrix* yT, * temp, *temp2,* extentedMatrix;
	float* lr;

	
	//[fliplr(x(:,1:ss)) x  x(:,n  :-1: n-p-s1+1)]
	temp = MatrixCut(x->mat, x->height, x->width, 0, x->height - 1, 0, ss - 1); // x(:,1:ss)
	
	temp2 = new Matrix(temp->height, temp->width);
	temp2->mat = Fliplr(temp->mat, temp->height, temp->width);// fliplr(x(:,1:ss))
	extentedMatrix = MatrixColExtend(temp2->mat, temp->height, temp->width, x->mat, x->height, x->width);
	delete temp; delete temp2;

	temp = MatrixCut(x->mat, x->height, x->width, 0, x->height - 1, n - 1, n - p - s1 + 1 - 1, 1, -1); // x(:, n : -1 : n - p - s1 + 1)
	yT = MatrixColExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);
	delete temp;

	//[flipud(yT(1:rr, : )); yT;  yT(m  :-1 : m - q - s2 + 1, : )]
	temp = MatrixCut(yT->mat, yT->height, yT->width, 0, rr - 1, 0, yT->width - 1); //yT(1:rr, : )
	temp2 = new Matrix(temp->height, temp->width);
	temp2->mat = Flipud(temp->mat, temp->height, temp->width);		//flipud(yT(1:rr, : ))
	delete temp; delete extentedMatrix;

	extentedMatrix = MatrixRowExtend(temp2->mat, temp2->height, temp2->width, yT->mat, yT->height, yT->width);
	delete temp2;
	
	temp = MatrixCut(yT->mat, yT->height, yT->width, m - 1, m - q - s2 + 1 - 1, 0, yT->width - 1, -1, 1);	//yT(m  :-1 : m - q - s2 + 1, : )
	
	delete yT;
	yT = MatrixRowExtend(extentedMatrix->mat, extentedMatrix->height, extentedMatrix->width, temp->mat, temp->height, temp->width);
	delete extentedMatrix; delete temp;

	// yT(1:m+p-1 ,1:n+q-1) 
	Matrix* retyT = MatrixCut(yT->mat, yT->height, yT->width, 0, m + p - 1 - 1, 0, n + q - 1 - 1);
	delete yT;

	return retyT;
}


/// <summary>
///		This function does not actually upsample the filter, 
///		it computes the convolution as if the filter had been upsampled.
///		This is the ultimate optimized version. Further optimized for separable(diagonal) upsampling matrices.
/// </summary>
/// <param name="signal : ">
///		A 2D Signal
/// </param>
/// <param name="filter : ">
///		2D filter
/// </param>
/// <param name="upMatrix : ">
///		separable upsampling matrix
/// </param>
/// <returns>
///		2D result of convolution with filter upsampled by a m, only the 'valid' part is returned. <para></para>
///		Similar to conv2(x, h, 'valid'), where h is the upsampled filter.
/// </returns>
Matrix* Atrousc(Matrix* signal, Matrix* filter, float* upMatrix)
{
	/* FArray   - Filter coefficients
	   SArray   - Signal coefficients
	   outArray - Output coefficients
	   M        - upsampling matrix 	*/
	float* FArray, * SArray, * outArray, * M;
	int SColLength, SRowLength, FColLength, FRowLength, O_SColLength, O_SRowLength;
	int SFColLength, SFRowLength;
	int n1, n2, k1, k2, f1, f2, kk2, kk1;
	float sum;
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
					f1 = SFRowLength - k1;	/* flipped index */
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


/// <summary>
///		ATROUSFILTERS	Generate pyramid 2D filters
/// </summary>
/// <param name="fname : ">
///		'maxflat':		Filters derived from 1-D using maximally flat mapping function with 4 vanishing moments
/// </param>
/// <returns>		 
///		h0, h1, g0, g1:	pyramid filters for 2-D nonsubsampled filter bank (lowpass and highpass)
/// </returns>
Cont* AtrousFilters(const char* fname) {


	int h0Width = 7;
	int h0Height = 7;

	int h1Width = 10;
	int h1Height = 10;

	int g0Width = 10;
	int g0Height = 10;

	int g1Width = 7;
	int g1Height = 7;

	if (!strcmp(fname, "maxflat")) {

		float h0[] = { -7.900496718847182e-07, 0., 0.000014220894093924927, 0.000025281589500310983, -0.000049773129328737247, -0.00022753430550279883, -0.00033182086219158167,
		   0,               0,              0,                   0,                       0,                          0,                     0,
		   0.000014220894093924927, 0., -0.0002559760936906487, -0.00045506861100559767,   0.0008959163279172705,   0.004095617499050379,    0.00597277551944847,
		   0.000025281589500310983, 0., -0.00045506861100559767, 0.0009765625             ,0.0015927401385195919, -0.0087890625, -0.01795090623402861,
		  -0.000049773129328737247, 0.,  0.0008959163279172705,  0.0015927401385195919, -0.0031357071477104465, -0.014334661246676327, -0.020904714318069645,
		  -0.00022753430550279883,  0.,  0.004095617499050379 ,-0.0087890625, -0.014334661246676327,    0.0791015625      ,      0.16155815610625748,
		  -0.00033182086219158167,  0.,  0.00597277551944847, -0.01795090623402861, -0.020904714318069645,    0.16155815610625748,     0.3177420190660832 };

		float g0[] = { -6.391587676622346e-010,             0.,                1.7257286726880333e-08,    3.067962084778726e-08, -1.3805829381504267e-07, -5.522331752601707e-07,
			-3.3747582932565985e-07,    1.9328161134105974e-06,     5.6949046198705095e-06,    7.649452131381623e-06,
			0.0,                          0.,                      0.,                          0.,               0.,  0.,  0.,  0.,  0.,  0.,
			 1.7257286726880333e-08,             0., -4.65946741625769e-07, -8.283497628902559e-07,    3.727573933006152e-06,    0.000014910295732024608,
			 9.111847391792816e-06, -0.000052186035062086126, -0.00015376242473650378, -0.00020653520754730382,
			 3.067962084778726e-08,              0., -8.283497628902559e-07, -1.2809236054493144e-06,   6.6267981031220475e-06,   0.00002305662489808766,
			 0.000010064497559808503, -0.0000806981871433068, -0.00021814634152337594, -0.00028666046030363884,
			-1.3805829381504267e-07,             0.,                3.727573933006152e-06,     6.6267981031220475e-06, -0.000029820591464049215, -0.00011928236585619686,
			-0.00007289477913434253,    0.000417488280496689,       0.0012300993978920302,     0.0016522816603784306,
			-5.522331752601707e-07,              0. ,               0.000014910295732024608 ,  0.00002305662489808766 ,-0.00011928236585619686, -0.00041501924816557786,
			-0.00018116095607655303,    0.0014525673685795225 ,     0.0039266341474207675 ,    0.005159888285465499,
			-3.3747582932565985e-07 ,            0. ,               9.111847391792816e-06,     0.000010064497559808503 ,-0.00007289477913434253, -0.00018116095607655303,
			 0.001468581806076247,     0.0006340633462679356, -0.01181401175635013, -0.021745034491193898,
			 1.9328161134105974e-06,             0., -0.000052186035062086126, -0.0000806981871433068,    0.000417488280496689,    0.0014525673685795225,
			 0.0006340633462679356, -0.005083985790028328, -0.013743219515972684, -0.018059608999129246,
			 5.6949046198705095e-06,             0., -0.00015376242473650378, -0.00021814634152337594,   0.0012300993978920302,   0.0039266341474207675 ,
			-0.01181401175635013, -0.013743219515972684 ,      0.0826466923977296  ,      0.1638988884584603,
			 7.649452131381623e-06 ,             0., -0.00020653520754730382, -0.00028666046030363884,   0.0016522816603784306,   0.005159888285465499 ,
			-0.021745034491193898 ,-0.018059608999129246  ,     0.1638988884584603 ,       0.31358726209239235 };

		float g1[] = { -7.900496718847182e-07,    0.,      0.000014220894093924927, 0.000025281589500310983, -0.000049773129328737247, -0.00022753430550279883, -0.00033182086219158167,
				0,                 0   ,             0      ,          0        ,                    0      ,                 0        ,               0,
			0.000014220894093924927,  0., -0.0002559760936906487, -0.00045506861100559767,   0.0008959163279172705,    0.004095617499050379 ,   0.00597277551944847,
			0.000025281589500310983,  0., -0.00045506861100559767, -0.0009765625 ,            0.0015927401385195919  ,  0.0087890625    ,        0.01329909376597139,
			-0.000049773129328737247,  0.,      0.0008959163279172705 ,  0.0015927401385195919, -0.0031357071477104465, -0.014334661246676327, -0.020904714318069645,
			-0.00022753430550279883,   0.,     0.004095617499050379,    0.0087890625, -0.014334661246676327, -0.0791015625, -0.1196918438937425,
			-0.00033182086219158167,   0.,      0.00597277551944847,     0.01329909376597139, -0.020904714318069645, -0.1196918438937425,      0.8177420190660831 };

		float h1[] = { 6.391587676622346e-010,       0.0, -1.7257286726880333e-08 ,-3.067962084778726e-08     ,1.3805829381504267e-07  ,5.522331752601707e-07 ,
			3.3747582932565985e-07, -1.9328161134105974e-06 ,-5.6949046198705095e-06, -7.649452131381623e-06,
			0.0,                     0.0,                           0.0,                            0.0,                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			-1.7257286726880333e-08       ,0.0,                  4.65946741625769e-07             ,8.283497628902559e-07 ,-3.727573933006152e-06 ,-0.000014910295732024608,
			-9.111847391792816e-06,   0.000052186035062086126  ,0.00015376242473650378           ,0.00020653520754730382,
			-3.067962084778726e-08        ,0.0,                  8.283497628902559e-07 ,-2.9917573832012203e-07 ,-6.6267981031220475e-06   ,5.3851632897621965e-06 ,
			0.00004049868144081346 ,-0.00001884807151416769 ,-0.00023692226948222173 ,-0.0003769812640795245,
			1.3805829381504267e-07       ,0.0, -3.727573933006152e-06 ,-6.6267981031220475e-06    ,0.000029820591464049215  ,0.00011928236585619686 ,
			0.00007289477913434253 ,-0.000417488280496689 ,-0.0012300993978920302 ,-0.0016522816603784306,
			5.522331752601707e-07        ,0.0, -0.000014910295732024608          ,5.3851632897621965e-06    ,0.00011928236585619686 ,-0.00009693293921571956 ,
			-0.0007289762659346422   ,0.00033926528725501844   ,0.004264600850679991             ,0.006785662753431441,
			3.3747582932565985e-07       ,0.0, -9.111847391792816e-06            ,0.00004049868144081346    ,0.00007289477913434253 ,-0.0007289762659346422 ,
			-0.001468581806076247    ,0.002551416930771248     ,0.01181401175635013              ,0.017093222023136675,
			-1.9328161134105974e-06       ,0.0,                  0.000052186035062086126 ,-0.00001884807151416769 ,-0.000417488280496689     ,0.00033926528725501844 ,
			0.002551416930771248 ,-0.0011874285053925643 ,-0.01492610297737997 ,-0.023749819637010044,
			-5.6949046198705095e-06       ,0.0,                  0.00015376242473650378 ,-0.00023692226948222173 ,-0.0012300993978920302    ,0.004264600850679991 ,
			0.01181401175635013 ,-0.01492610297737997 ,-0.0826466923977296 ,-0.12203257624594532,
			-7.649452131381623e-06        ,0.0,                  0.00020653520754730382 ,-0.0003769812640795245 ,-0.0016522816603784306    ,0.006785662753431441 ,
			0.017093222023136675 ,-0.023749819637010044 ,-0.12203257624594532               ,0.821896776039774 };


		/*
			g0 = [g0   fliplr(g0(:, 1 : end - 1))];
			g0 = [g0; flipud(g0(1:end - 1, : ))];
			h0 = [h0   fliplr(h0(:, 1 : end - 1))];
			h0 = [h0; flipud(h0(1:end - 1, : ))];

			g1 = [g1   fliplr(g1(:, 1 : end - 1))];
			g1 = [g1; flipud(g1(1:end - 1, : ))];
			h1 = [h1   fliplr(h1(:, 1 : end - 1))];
			h1 = [h1; flipud(h1(1:end - 1, : ))];
		*/

		float* lr, *ud;
		// g0
		Matrix* cutG0 = MatrixCut(g0, g0Height, g0Width, 0, g0Height - 1, 0, g0Width - 2);
		lr = Fliplr(cutG0->mat, cutG0->height, cutG0->width);
		Matrix* extG0 = MatrixColExtend(g0, g0Height, g0Width, lr, cutG0->height, cutG0->width);
		delete[] lr;

		// g0
		Matrix* cutExtG0 = MatrixCut(extG0->mat, extG0->height, extG0->width, 0, extG0->height - 2, 0, extG0->width - 1);
		ud = Flipud(cutExtG0->mat, cutExtG0->height, cutExtG0->width);
		Matrix* retG0 = MatrixRowExtend(extG0->mat, extG0->height, extG0->width, ud , cutExtG0->height, cutExtG0->width);
		delete[] ud;

		// h0
		Matrix* cutH0 = MatrixCut(h0, h0Height, h0Width, 0, h0Height - 1, 0, h0Width - 2);
		lr = Fliplr(cutH0->mat, cutH0->height, cutH0->width);
		Matrix* extH0 = MatrixColExtend(h0, h0Height, h0Width, lr, cutH0->height, cutH0->width);
		delete[] lr; 

		// h0
		Matrix* cutExtH0 = MatrixCut(extH0->mat, extH0->height, extH0->width, 0, extH0->height - 2, 0, extH0->width - 1);
		ud = Flipud(cutExtH0->mat, cutExtH0->height, cutExtH0->width);
		Matrix* retH0 = MatrixRowExtend(extH0->mat, extH0->height, extH0->width, ud , cutExtH0->height, cutExtH0->width);
		delete[] ud;

		// g1
		Matrix* cutG1 = MatrixCut(g1, g1Height, g1Width, 0, g1Height - 1, 0, g1Width - 2);
		lr = Fliplr(cutG1->mat, cutG1->height, cutG1->width);
		Matrix* extG1 = MatrixColExtend(g1, g1Height, g1Width, lr, cutG1->height, cutG1->width);
		delete[] lr;

		// g1
		Matrix* cutExtG1 = MatrixCut(extG1->mat, extG1->height, extG1->width, 0, extG1->height - 2, 0, extG1->width - 1);
		ud = Flipud(cutExtG1->mat, cutExtG1->height, cutExtG1->width);
		Matrix* retG1 = MatrixRowExtend(extG1->mat, extG1->height, extG1->width, ud , cutExtG1->height, cutExtG1->width);
		delete[] ud;

		// h1
		Matrix* cutH1 = MatrixCut(h1, h1Height, h1Width, 0, h1Height - 1, 0, h1Width - 2);
		lr = Fliplr(cutH1->mat, cutH1->height, cutH1->width);
		Matrix* extH1 = MatrixColExtend(h1, h1Height, h1Width, lr, cutH1->height, cutH1->width);
		delete[] lr;

		// h1
		Matrix* cutExtH1 = MatrixCut(extH1->mat, extH1->height, extH1->width, 0, extH1->height - 2, 0, extH1->width - 1);
		ud = Flipud(cutExtH1->mat, cutExtH1->height, cutExtH1->width);
		Matrix* retH1 = MatrixRowExtend(extH1->mat, extH1->height, extH1->width, ud , cutExtH1->height, cutExtH1->width);
		delete[] ud;


		Cont* ret = new Cont(4);
		ret->mats[0] = retG0;
		ret->mats[1] = retH0;
		ret->mats[2] = retG1;
		ret->mats[3] = retH1;

		delete cutG0; delete extG0; delete cutExtG0;
		delete cutH0; delete extH0; delete cutExtH0;
		delete cutG1; delete extG1; delete cutExtG1;
		delete cutH1; delete extH1; delete cutExtH1;

		return ret;
	}

	else exit(1);
}


/// <summary>
///		ATROUSDEC - computes the 2 - D atrous decomposition using symmetric extension.
/// </summary>
/// <param name="image : ">
///		input image
/// </param>
/// <param name="filters : ">
///		can be any filter available in the function atrousfilters
/// </param>
/// <param name="level : ">
///		N levels - number of decomposition levels
/// </param>
/// <returns>
///		y->[0] => LFC (Low frequency coefficients 
///		<para></para>
///		y->[1..4] => HFC (High frequency coefficients)
/// </returns>
Cont* AtrousDec(Matrix* image, Cont* filters, int level)
{
	Matrix* y0, *y1,* sym, * upSample;

	float* shift = new float[2]{ 1.0F, 1.0F };

	// NSLP - Non Subsampled Laplacian Pyramid
	// Orijinal goruntuye alcak ve yuksek katsayilarin uygulanmasi 
	// symext kullanilarak katsayilara gore genisletilir, ardindan conv2 ile kucultulurek ayni boyuta getirilir.
	sym = symext(image, filters->mats[1], shift);
	y0 = Conv2(sym, filters->mats[1], "valid");	// Alcak
	delete sym;

	sym = symext(image, filters->mats[3], shift);
	y1 = Conv2(sym, filters->mats[3], "valid");	// Yuksek
	delete sym;

	Cont* y = new Cont(level + 1);
	y->mats[level] = y1;
	image = y0;

	float* I2 = Eye(2), * I2L;
	int L;

	// Seviye 1'e inene kadar her alcak frekans katsayilarinin alcak ve yuksek katsayyisi hesaplanir.
	// Yuksek frekansa katsayilari y[4..1] olarak tutulur.
	for (int i = 1; i <= level - 1; i++)
	{
		shift[0] = (-1.0F * (float)pow(2, i - 1)) + 2.0F;
		shift[1] = (-1.0F * (float)pow(2, i - 1)) + 2.0F;
		L = (int)pow(2, i);
		I2L = ScalarMatMul(I2, 2 * 2, L);
		
		upSample = Upsample2df(filters->mats[1], i);
		sym = symext(image, upSample, shift);
		y0 = Atrousc(sym, filters->mats[1], I2L);
		delete upSample;
		delete sym;

		upSample = Upsample2df(filters->mats[3], i);
		sym = symext(image, upSample, shift);
		y1 = Atrousc(sym, filters->mats[3], I2L);
		delete upSample;
		delete sym;

		y->mats[level - i] = y1;
		image = y0;

		delete[] I2L; 
	}
	y->mats[0] = image;		// Alcak frekans katsayisi atanir.

	delete[] I2;
	delete[] shift;

	return y;
}


/// <summary>
///		SATROUSREC - computes the inverse of 2 - D atrous decomposition computed with ATROUSDEC
///		y = satrousdrc(x, fname)
/// </summary>
/// <param name="y : ">
///		image
/// </param>
/// <param name="lpfilt : ">
///		 can be any filter available in the function atrousfilters
///	</param>
/// <returns>
///		x : reconstructed image
/// </returns>
Matrix* AtrousRec(Cont* y, const char* lpfilt) {

	int NLevels = y->matNums - 1;

	// [g0, h0, g1, h1] = atrousfilters(fname);
	Cont* ret = AtrousFilters(lpfilt);
	Matrix* g0 = ret->mats[0];
	Matrix* h0 = ret->mats[1];
	Matrix* g1 = ret->mats[2];
	Matrix* h1 = ret->mats[3];


	Matrix* x;
	Matrix* y1;

	x = y->mats[0];

	Matrix* I2 = EyeMatrix(2);
	float* shift = new float[2]{ 1.0, 1.0 };
	float L = 0.0;

	Matrix* up, * sym, * mult, * atrousc1, * atrousc2;
	for (int i = NLevels - 1; i >= 1; i--) {

		y1 = y->mats[NLevels - i];			// Matlab: y{2} <=> C: y[1]

		shift[0] = -1 * pow(2, (i - 1)) + 2.0;
		shift[1] = -1 * pow(2, (i - 1)) + 2.0;

		L = pow(2, i);


		// x = atrousc(symext(x, upsample2df(g0, i), shift), g0, L * I2) + atrousc(symext(y1, upsample2df(g1, i), shift), g1, L * I2);

		mult = *I2 * L;

		up = Upsample2df(g0, i);
		sym = symext(x, up, shift);		delete up;
		atrousc1 = Atrousc(sym, g0, mult->mat);	delete sym;

		up = Upsample2df(g1, i);
		sym = symext(y1, up, shift);	delete up;
		atrousc2 = Atrousc(sym, g1, mult->mat);	delete sym;

		x = *atrousc1 + *atrousc2;


		 delete mult; delete atrousc1; delete atrousc2;
	}

	shift[0] = 1.0;
	shift[1] = 1.0;

	Matrix* symetxX = symext(x, g0, shift);
	Matrix* symetxY = symext(y->mats[NLevels], g1, shift);

	Matrix* conv2X = Conv2(symetxX, g0, "valid");
	Matrix* conv2Y = Conv2(symetxY, g1, "valid");
	x = *conv2X + *conv2Y;

	delete I2; delete g0; delete h0; delete g1; delete h1; delete symetxX; delete symetxY; delete conv2X; delete conv2Y;

	return x;
}