#pragma once
#include <stdlib.h>

class Matrix {

public:
	int height;
	int width;
	int depth;
	float* mat;

	Matrix(){}
	~Matrix() { delete[] mat; }
	
	Matrix(int height, int width, int depth = 1) {
		this->height = height;
		this->width = width;
		this->depth = depth;
	}

	void CreateMatrix(int height, int width, int depth = 1) {

		this->height = height;
		this->width = width;
		this->depth = depth;
		mat = new float[height * width * depth]();
	}

	int GetSize2D() { return height * width; }
	int GetSize3D() { return height * width * depth; }


	Matrix* operator+(Matrix& b) {

		Matrix* temp = new Matrix;
		temp->CreateMatrix(b.height, b.width, b.depth);

		for (int i = 0; i < temp->GetSize3D(); i++)
		{
			temp->mat[i] = this->mat[i] + b.mat[i];
		}

		return temp;
	}

	Matrix* operator*(float L) {
		
		Matrix* temp = new Matrix;
		temp->CreateMatrix(this->height, this->width, this->depth);

		for (int i = 0; i < temp->GetSize3D(); i++)
		{
			temp->mat[i] = this->mat[i] * L;
		}
		
		return temp;
	}

	void operator/=(float L) {

		for (int i = 0; i < this->GetSize3D(); i++)
		{
			this->mat[i] = this->mat[i] / L;
		}
	}

};

// Cont => Container class
class Cont {

public:
	Matrix** mats;
	int matNums;

	Cont() {}

	Cont(int cellNums) {
		this->matNums = cellNums;
		mats = new Matrix*[cellNums];
	}

	void CreateCells(int index, int size) {
		mats[index] = new Matrix[size];	
	}

	~Cont(){

		int depth;
		for (size_t cell = 0; cell < matNums; cell++) {

			depth = mats[cell]->depth;
			for (size_t d = 0; d < depth; d++)
				delete[] mats[cell][d].mat;
		}
		delete[] mats;
	}
};