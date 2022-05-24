//
// Created by fatihbasatemur on 4/30/2022.
//

#ifndef NSST_FFTWF_WRAPPER_H
#define NSST_FFTWF_WRAPPER_H

#include <Eigen/Core>
#include "fftwf_.h"
#include <tuple>
#include <map>

class ForwEigen : public ForwBack<ArrayXXf> {
public:
    void SetIn(const ArrayXXf *real, const ArrayXXf *imag) final {
        // ColMajor to RowMajor settle
        int in_index = 0;

        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                in_index = r * cols + c;
                in[in_index][REAL] = real->coeff(r, c);
                in[in_index][IMAG] = imag->coeff(r, c);
            }
    }

    void GetOut(ArrayXXf *real, ArrayXXf *imag) final {
        // RowMajor to ColMajor getter
        float *real_p = real->data();
        float *imag_p = imag->data();
        int _size = rows * cols;

        for (int i = 0; i < _size; i++) {
            real_p[i] = out[i][REAL];
            imag_p[i] = out[i][IMAG];
        }
    }

    friend class Eigen_2D;
};

class BackEigen : public ForwBack<ArrayXXf> {
protected:
    float normSize = 0.0F;
public:
    void SetIn(const ArrayXXf *real, const ArrayXXf *imag) final {
        // ColMajor to RowMajor settle
        int in_index = 0;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                in_index = r * cols + c;
                in[in_index][REAL] = real->coeff(r, c);
                in[in_index][IMAG] = imag->coeff(r, c);
            }
    }

    void SetIn(const ArrayXXf &real) {
        // ColMajor to RowMajor settle
        int in_index = 0;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                in_index = r * cols + c;
                in[in_index][REAL] = real.coeff(r, c);
                in[in_index][IMAG] = 0;
            }
    }

    void GetOut(ArrayXXf* real, ArrayXXf* imag) final {
        // RowMajor to ColMajor getter

        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                real[0](r,c) = out[r * rows + c ][REAL] / normSize;
                imag[0](r,c) = out[r * rows + c ][IMAG] / normSize;
            }
    }
    // TODO: getterlarda non-contiguous array mapping kullanilabilir
    void GetOut(ArrayXXf &real) {

        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                real(r,c) = out[r * rows + c ][REAL] / normSize;
            }
    }

    friend class Eigen_2D;
};

class Eigen_2D : public FFTWF<ForwEigen, BackEigen> {
public:
private:
    void CreatePlan() override {
        forward->plan = fftwf_plan_dft_2d(forward->rows, forward->cols, forward->in, forward->out, FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        backward->plan = fftwf_plan_dft_2d(backward->rows, backward->cols, backward->in, backward->out, FFTW_BACKWARD,
                                           FFTW_ESTIMATE);
        backward->normSize = static_cast<float>(backward->size);
    }

public:
    Eigen_2D(int rows, int cols) {
        forward->rows = backward->rows = rows, forward->cols = backward->cols = cols, forward->size = backward->size =
                rows * cols;
        Alloc();
        Eigen_2D::CreatePlan();
    }
};

class FFTWF_Wrapper{
private:
    // use Flyweight pattern
    inline static std::map <std::tuple<int, int>, FFTWF<ForwEigen, BackEigen>*> _store_map;
    inline static FFTWF<ForwEigen, BackEigen>* _iter = nullptr;

public:
    class Eigen{
    public:
        class _2D{
        public:
            static FFTWF<ForwEigen, BackEigen>* Get(int rows, int cols) {

                auto search = _store_map.find(std::tuple(rows, cols));
                if (search != _store_map.end()) {
                    _iter = search->second;
                }
                else {
                    _iter = new Eigen_2D(rows, cols);
                    _store_map[std::tuple(rows, cols)] = _iter;
                }
                return _iter;
            }

            ~_2D(){
                for(auto it = _store_map.begin(); it != _store_map.end(); ++it)
                    delete it->second;
                _store_map.clear();
                _iter = nullptr;
            }
        };
    };
};

#endif //NSST_FFTWF_WRAPPER_H
