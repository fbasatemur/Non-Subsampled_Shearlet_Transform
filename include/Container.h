#ifndef NSST_CONTAINER_H
#define NSST_CONTAINER_H

#include <Eigen/Core>
#include "Definitions.h"

class Tensor {

    class Create{
    public:
        inline static Tensor* _self = nullptr;
        static void Default(){
            _self->_mat = ArrayXXf(_self->_h * _self->_d, _self->_w);
        }
        static void Zero(){
            _self->_mat = ArrayXXf::Zero(_self->_h * _self->_d, _self->_w);
        }
    };

public:
    int _w = 0, _h = 0, _d = 1;
    ArrayXXf _mat;

    Tensor(){}

    Tensor(Tensor* in){
        in->CopyTo(this);
    }

    Create& Set(int h, int w, int d = 1){
        _h = h, _w = w, _d = d;
        Create::_self = this;
        static Create _create;
        return _create;
    }

    void CopyTo(Tensor* tensor){
        tensor->Set(_h, _w, _d);
        if (nullptr == tensor->_mat.data())
            Create::Default();
        tensor->_mat = this->_mat;
    }

    Tensor& operator+(const Tensor& in) {
        Tensor* ret = new Tensor;
        ret->Set(this->_h, this->_w, this->_d).Default();
        ret->_mat = this->_mat + in._mat;
		return *ret;
	}

    Tensor& operator/(const Tensor& in) {
        Tensor* ret = new Tensor;
        ret->Set(this->_h, this->_w, this->_d).Default();
        ret->_mat = this->_mat / in._mat;
        return *ret;
    }

    Tensor& operator*(float L) {
        Tensor* ret = new Tensor;
        ret->Set(this->_h, this->_w, this->_d).Default();
        ret->_mat = this->_mat * L;
		return *ret;
	}

	void operator/=(float L) {
        this->_mat = (this->_mat / L).eval();
	}

    void operator=(ArrayXXf _mat){
        this->Set(_mat.rows(), _mat.cols());
        this->_mat = _mat;
    }

};

// Cont => Container class
class Cont {

public:
	Tensor** _conts;
	int _cont_num;

	Cont(int cont_num) {
		this->_cont_num = cont_num;
		_conts = new Tensor*[cont_num];
	}

    Tensor* Create_Cell(int cont_index, int depth) {
        _conts[cont_index] = new Tensor[depth];
        return _conts[cont_index];
    }

	~Cont(){
        for (int i = 0; i < _cont_num; ++i) {
            _conts[i]->_d == 1 ? delete _conts[i] : delete[] _conts[i];
        }
		delete[] _conts;
	}
};

#endif