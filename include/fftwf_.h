//
// Created by fbasatemur on 5/1/2022.
//

#ifndef NSST_FFTWF_H
#define NSST_FFTWF_H

#include "fftw-3.3.5-dll64/fftw3.h"
#define REAL 0
#define IMAG 1

class BackComplex;
template <typename T_Forward, typename T_Backward> class FFTWF;

template<typename T>
class ForwBack {
protected:
    int rows = 1, cols = 1, depth = 1, size = 1;
    fftwf_plan plan = nullptr;
    fftwf_complex *in = nullptr, *out = nullptr;

public:
    virtual void SetIn(const T *, const T *) = 0;
    virtual void GetOut(T *, T *) = 0;

    void Alloc() {
        in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * size);
        out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * size);
    }

    void Execute() { fftwf_execute(plan); };

    virtual ~ForwBack() {
        fftwf_destroy_plan(plan);
        fftwf_free(in);
        fftwf_free(out);
    }
};

class ForwComplex : public ForwBack<float> {
public:
    void SetIn(const float * real, const float * imag) final {
        for (int i = 0; i < size; i++) {
            in[i][REAL] = real[i];
            in[i][IMAG] = imag[i];
        }
    }

    void GetOut(float * real, float * imag) final {
        for (int i = 0; i < size; i++) {
            real[i] = out[i][REAL];
            imag[i] = out[i][IMAG];
        }
    }

    friend class FFTWF<ForwComplex, BackComplex>;
    friend class Complex_1D;
    friend class Complex_2D;
    friend class Complex_2D_Each_Row;
    friend class Complex_2D_Each_Col;
    friend class Complex_3D_Each_Channel;
};

class BackComplex : public ForwBack<float> {
protected:
    int normSize = 0;

public:
    void SetIn(const float * real, const float * imag) final {
        for (int i = 0; i < size; i++) {
            in[i][REAL] = real[i];
            in[i][IMAG] = imag[i];
        }
    }

    void GetOut(float * real, float * imag) final {
        for (int i = 0; i < size; i++) {
            real[i] = out[i][REAL] / float(normSize);
            imag[i] = out[i][IMAG] / float(normSize);
        }
    }

    friend class FFTWF<ForwComplex, BackComplex>;
    friend class Complex_1D;
    friend class Complex_2D;
    friend class Complex_2D_Each_Row;
    friend class Complex_2D_Each_Col;
    friend class Complex_3D_Each_Channel;
};

template<typename T_Forward = ForwComplex, typename T_Backward = BackComplex>
class FFTWF {
protected:
    T_Forward *forward = new T_Forward;
    T_Backward *backward = new T_Backward;

    virtual void CreatePlan() = 0;

    void Alloc() {
        forward->Alloc();
        backward->Alloc();
    }

public:
    T_Forward *Forward() { return forward; }
    T_Backward *Backward() { return backward; }

    ~FFTWF() {
        delete forward, delete backward;
    }
};

class Complex_1D : public FFTWF<ForwComplex, BackComplex> {
private:
    void CreatePlan() override {
        forward->plan = fftwf_plan_dft_1d(forward->size, forward->in, forward->out, FFTW_FORWARD, FFTW_ESTIMATE);
        backward->plan = fftwf_plan_dft_1d(backward->size, backward->in, backward->out, FFTW_BACKWARD, FFTW_ESTIMATE);
        backward->normSize = backward->size;
    }

public:
    Complex_1D(int size) {
        forward->cols = backward->cols = size, forward->size = backward->size = size;
        Alloc();
        CreatePlan();
    }
};

class Complex_2D : public FFTWF<ForwComplex, BackComplex> {
private:
    void CreatePlan() override {
        forward->plan = fftwf_plan_dft_2d(forward->rows, forward->cols, forward->in, forward->out, FFTW_FORWARD, FFTW_ESTIMATE);
        backward->plan = fftwf_plan_dft_2d(backward->rows, backward->cols, backward->in, backward->out, FFTW_BACKWARD, FFTW_ESTIMATE);
        backward->normSize = backward->size;
    }

public:
    Complex_2D(int rows, int cols) {
        forward->rows = backward->rows = rows, forward->cols = backward->cols = cols, forward->size = backward->size = rows * cols;
        Alloc();
        CreatePlan();
    }
};

class Complex_2D_Each_Row : public FFTWF<ForwComplex, BackComplex> {
private:
    void CreatePlan() override {
        /*
            in[] -> *(in + col_index * istride + height * idist)

            if data is column-major (fft for each column), set istride=howmany, idist=1
            if data is row-major	(fft for each row), set istride=1, idist=N
        */
        int rank = 1;					// we are computing 1d transforms
        int n[] = { forward->cols };	// 1d transforms of length
        int howmany = forward->rows, * inembed = n, * onembed = n, idist, odist, istride, ostride;
        idist = odist = forward->cols;
        istride = ostride = 1;			// distance between two elements in the same column
        forward->plan = fftwf_plan_many_dft(rank, n, howmany, forward->in, inembed, istride, idist, forward->out, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        backward->plan = fftwf_plan_many_dft(rank, n, howmany, backward->in, inembed, istride, idist, backward->out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
        backward->normSize = backward->cols;
    }

public:
    Complex_2D_Each_Row(int rows, int cols) {
        forward->rows = backward->rows = rows, forward->cols = backward->cols = cols, forward->size = backward->size = rows * cols;
        Alloc();
        CreatePlan();
    }
};

class Complex_2D_Each_Col : public FFTWF<ForwComplex, BackComplex> {
private:
    void CreatePlan() override {
        /*
            in[] -> *(in + col_index * istride + height * idist)

            if data is column-major (fft for each column), set istride=howmany, idist=1
            if data is row-major	(fft for each row), set istride=1, idist=N
        */
        int rank = 1;					// we are computing 1d transforms
        int n[] = { forward->rows };	// 1d transforms of length
        int howmany = forward->cols, * inembed = n, * onembed = n, idist, odist, istride, ostride;
        idist = odist = 1;
        istride = ostride = forward->cols;			// distance between two elements in the same column
        forward->plan = fftwf_plan_many_dft(rank, n, howmany, forward->in, inembed, istride, idist, forward->out, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        backward->plan = fftwf_plan_many_dft(rank, n, howmany, backward->in, inembed, istride, idist, backward->out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
        backward->normSize = backward->rows;
    }

public:
    Complex_2D_Each_Col(int rows, int cols) {
        forward->rows = backward->rows = rows, forward->cols = backward->cols = cols, forward->size = backward->size = rows * cols;
        Alloc();
        CreatePlan();
    }

};

class Complex_3D_Each_Channel : public FFTWF<ForwComplex, BackComplex> {
private:
    void CreatePlan() override {
        int rank = 2;					// we are computing 2d transforms
        int n[] = { forward->rows, forward->cols };	// 2d transforms of size
        int howmany = forward->depth, * inembed = n, * onembed = n, idist, odist, istride, ostride;
        idist = odist = forward->rows * forward->cols;
        istride = ostride = 1;
        forward->plan = fftwf_plan_many_dft(rank, n, howmany, forward->in, inembed, istride, idist, forward->out, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        backward->plan = fftwf_plan_many_dft(rank, n, howmany, backward->in, inembed, istride, idist, backward->out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
        backward->normSize = backward->rows * backward->cols;
    }

public:
    Complex_3D_Each_Channel(int rows, int cols, int depth) {
        forward->rows = backward->rows = rows, forward->cols = backward->cols = cols, forward->depth = backward->depth = depth, forward->size = backward->size = rows * cols * depth;
        Alloc();
        CreatePlan();
    }
};


class FFTWF_Factory {
public:
    class Complex {
    public:
        class _1D {
        public:
            static FFTWF<ForwComplex, BackComplex>* Create(int size) {
                return new Complex_1D(size);
            }
        };
        class _2D {
        public:
            static FFTWF<ForwComplex, BackComplex>* Create(int rows, int cols) {
                return new Complex_2D(rows, cols);
            }
        };
        class _2D_Each_Row {
        public:
            static FFTWF<ForwComplex, BackComplex>* Create(int rows, int cols) {
                return new Complex_2D_Each_Row(rows, cols);
            }
        };
        class _2D_Each_Col {
        public:
            static FFTWF<ForwComplex, BackComplex>* Create(int rows, int cols) {
                return new Complex_2D_Each_Col(rows, cols);
            }
        };
        class _3D_Each_Channel {
        public:
            static FFTWF<ForwComplex, BackComplex>* Create(int rows, int cols, int depth) {
                return new Complex_3D_Each_Channel(rows, cols, depth);
            }
        };
    };
};

#endif