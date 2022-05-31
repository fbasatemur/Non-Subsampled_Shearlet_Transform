#ifndef NSST_NSST_H
#define NSST_NSST_H

#include "Container.h"
#include "../Config.h"
#include <iostream>
#define PI 3.14159265358979323846F
#define max(a, b) (a > b ? a : b)

class NSST {
public:
    enum conv_type{full, same, valid};

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
    Tensor *GenXYCoordinates(int n);

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
    void ShearingFiltersMyer();

    /// <summary>
    ///		ATROUSFILTERS	Generate pyramid 2D filters
    ///		'maxflat':		Filters derived from 1-D using maximally flat mapping function with 4 vanishing moments
    ///		h0, h1, g0, g1:	pyramid filters for 2-D nonsubsampled filter bank (lowpass and highpass)
    /// </summary>
    void AtrousFilters();

    /// <summary>
    ///		This function generates the matrix that contains the number of times the polar grid points go through the rectangular grid
    /// </summary>
    /// <param name="L : ">
    ///		L is the order of the block matrix
    /// </param>
    /// <returns>
    ///		D is the common grid point values
    /// </returns>
    Tensor *AvgPol(int L, const ArrayXXf &x1, const ArrayXXf &y1, const ArrayXXf &x2, const ArrayXXf &y2);

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
    Tensor &RecFromPol(const ArrayXXf &l, int n, const Tensor* gen);

    /// <summary>
    ///		This function computes the(local) nonsubsampled shearlet transform as given.
    /// </summary>
    /// <param name="input">
    ///		input input
    /// </param>
    /// <returns>
    ///		decomposition of input
    /// </returns>
    Cont *Dec(Tensor *input);

    /// <summary>
    ///		ATROUSDEC - computes the 2 - D atrous decomposition using symmetric extension.
    /// </summary>
    /// <param name="input : ">
    ///		input input
    /// </param>
    /// <param name="level : ">
    ///		N levels - number of decomposition levels
    /// </param>
    /// <returns>
    ///		y->[0] => LFC (Low frequency coefficients)     <=> (AFK)
    ///		<para></para>
    ///		y->[1..4] => HFC (High frequency coefficients) <=> (YFK)
    /// </returns>
    Cont *AtrousDec(Tensor *input, int level);

    /// <summary>
    ///		Performs symmetric extension for input x, filter h. <para></para>
    ///		The filter h is assumed have odd dimensions. <para></para>
    ///		If the filter has horizontaland vertical symmetry, then the nonsymmetric part of conv2(h, x) has the same size of x.
    /// </summary>
    /// <param name="x : ">
    ///		mxn input
    /// </param>
    /// <param name="h : ">
    ///		2 - D filter coefficients
    /// </param>
    /// <param name="shift : ">
    ///		optional shift
    /// </param>
    /// <returns>
    ///		yT input symetrically extended(H / V symmetry)
    /// </returns>
    Tensor *Symext(const Tensor *x, const Tensor *h, const float *shift);

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
    Tensor *Atrousc(const Tensor *signal, const Tensor *filter, const Eigen::Array22f &upMatrix);

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
    ///		the reconstructed input
    /// </returns>
    Tensor *Rec(const Cont * __restrict dst);

    /// <summary>
    ///		SATROUSREC - computes the inverse of 2 - D atrous decomposition computed with ATROUSDEC
    ///		y = satrousdrc(x, fname)
    /// </summary>
    /// <param name="y : ">
    ///		input
    /// </param>
    /// <param name="lpfilt : ">
    ///		 can be any filter available in the function atrousfilters
    ///	</param>
    /// <returns>
    ///		x : reconstructed input
    /// </returns>
    Tensor *AtrousRec(const Cont *y);

    ~NSST() {
        for (int i = 0; i < Conf::sp.dcompSize; ++i) {
            delete[] shearing_filters[i];
        }
        delete[] shearing_filters;
    }

private:
    Tensor *Windowing(const float *x, int N, int L);

    Cont *atrous_filters = nullptr;
    Tensor **shearing_filters = nullptr;

    inline float *ONE(int N) {
        float *ret = new float[N];
        std::fill_n(ret, N, 1.0F);
        return ret;
    }

    inline float *ZERO(int N) {
        return new float[N]();
    }

    inline float Realmod(float x, float y) {
        float result = fmodf(x, y);
        return result >= 0 ? result : result + y;
    }

    inline float *Linspace(int d1, int d2, int N) {

        int n1 = N - 1;
        float a1 = d1, a2 = d2;
        float* y = ZERO(N);

        for (int i = 0; i <= n1; i++)
            y[i] = a1 + i * (a2 - a1) / n1;

        y[0] = a1; y[n1] = a2;

        return y;
    }

    inline Tensor *Upsample2df(const Tensor *h, int power) {
        int scale_k = pow(2, power);
        int height = scale_k * h->_h, width = scale_k * h->_w;

        Tensor *ho = new Tensor();
        ho->Set(height, width).Create::Zero();
        float* ho_p = ho->_mat.data();
        const float* h_p = h->_mat.data();

        int hStep = 0;
        for (int row = 0; row < height; row += scale_k) {
            for (int col = 0; col < width; col += scale_k) {
                ho_p[row * width + col] = h_p[hStep];
                hStep++;
            }
        }

        return ho;
    }

    // Input-Side convolution
    Tensor* Conv2D(const ArrayXXf& input, const ArrayXXf& kernel, conv_type type)
    {
        int outRow, outCol, edgeRows, edgeCols;
        int input_h = input.rows(), input_w = input.cols();
        int kernel_h = kernel.rows(), kernel_w = kernel.cols();

        switch (type) {
            case conv_type::full:
                outRow = input_h + kernel_h - 1;
                outCol = input_w + kernel_h - 1;
                edgeRows = kernel_h - 1;
                edgeCols = kernel_w - 1;
                break;

            case conv_type::same:
                outRow = input_h;
                outCol = input_w;
                edgeRows = (kernel_h - 1) / 2;
                edgeCols = (kernel_w - 1) / 2;
                break;

            case conv_type::valid:
                outRow = input_h - kernel_h + 1;
                outCol = input_w - kernel_w + 1;
                edgeRows = edgeCols = 0;
                break;

            default:
                exit(1);
        }

        Tensor* outMat = new Tensor;
        outMat->Set(outRow, outCol).Default();
        float* out_p = outMat->_mat.data();
        const float* __restrict input_p = input.data();
        const float* __restrict kernel_p = kernel.data();

        int iImage, iKernel, jImage, jKernel, jKernelTemp, jImageTemp;
        float sum;

        for (int i = 0; i < outRow; i++) {
            for (int j = 0; j < outCol; j++)
            {
                sum = 0.0F;
                iKernel = kernel_h - 1 - max(0, edgeRows - i);
                iImage = max(0, i - edgeRows);
                jKernelTemp = kernel_w - 1 - max(0, edgeCols - j);
                jImageTemp = max(0, j - edgeCols);

                for (; (iKernel >= 0) && (iImage < input_h); iKernel--, iImage++)
                    for (jKernel = jKernelTemp, jImage = jImageTemp; (jKernel >= 0) && (jImage < input_w); jKernel--, jImage++)
                        sum += input_p[input_w * iImage + jImage] * kernel_p[kernel_w * iKernel + jKernel];

                out_p[i * outCol + j] = sum;
            }
        }
        return outMat;
    }

    Tensor *Sum_Axis3(const Tensor *tensor) {

        Tensor *ret = new Tensor;
        ret->Set(tensor->_h, tensor->_w).Create::Zero();

        for (int d = 0; d < tensor->_d; d++)
            ret->_mat += tensor[d]._mat;

        return ret;
    }
};

#endif