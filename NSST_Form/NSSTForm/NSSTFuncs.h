#pragma once
#include "Container.h"

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
Matrix* GenXYCoordinates(int n);

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
Matrix* ShearingFiltersMyer(int n, int level);

// This function computes the Meyer window function of signal x
float	MeyerWind(float x);

// This function computes L band sequence of a Meyer windowed signal.
Matrix* Windowing(float* x, int lenghtX, int L);

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
Matrix* RecFromPol(Matrix* l, int n, Matrix* gen);

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
Matrix* Symext(Matrix* x, const Matrix* h, float* shift);

// upsample filter by 2^power
Matrix* Upsample2df(const Matrix* mat, int power);

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
Matrix* Atrousc(Matrix* signal, const Matrix* filter, float* upMatrix);

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
Cont* AtrousDec(Matrix* image, Cont* filters, int level);

/// <summary>
///		ATROUSFILTERS	Generate pyramid 2D filters
/// </summary>
/// <param name="fname : ">
///		'maxflat':		Filters derived from 1-D using maximally flat mapping function with 4 vanishing moments
/// </param>
/// <returns>		 
///		h0, h1, g0, g1:	pyramid filters for 2-D nonsubsampled filter bank (lowpass and highpass)
/// </returns>
Cont* AtrousFilters(const char* fname);

/// <summary>
///		SATROUSREC - computes the inverse of 2 - D atrous decomposition computed with ATROUSDEC
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
Matrix* AtrousRec(Cont* y, const Cont* filters);

/// <summary>
///		This function generates the matrix that contains the number of times the polar grid points go through the rectangular grid 
/// </summary>
/// <param name="L : ">
///		L is the order of the block matrix
/// </param>
/// <returns>
///		D is the common grid point values
/// </returns>
float* AvgPol(int L, float* x1, float* y1, float* x2, float* y2);