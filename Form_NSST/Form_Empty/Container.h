#pragma once
#include <stdlib.h>

class Matrix {

public:
	int height;
	int width;
	int depth;

	double* mat;

	Matrix(){}
	~Matrix(){}
	
	Matrix(int height, int width, int depth = 1) {
		this->height = height;
		this->width = width;
		this->depth = depth;
	}

	void CreateMatrix(int height, int width, int depth = 1) {

		this->height = height;
		this->width = width;
		this->depth = depth;
		mat = new double[height * width * depth];
	}

	int GetSize2D() { return height * width; }
	int GetSize3D() { return height * width * depth; }

	double operator()(int row, int col) {
		return this->mat[row * width + col];
	}

	double operator()(int index) {
		return this->mat[index];
	}
	
	Matrix* operator+(Matrix& b) {

		Matrix* temp = new Matrix;
		temp->CreateMatrix(b.height, b.width, b.depth);

		for (int i = 0; i < temp->GetSize3D(); i++)
		{
			temp->mat[i] = this->mat[i] + b.mat[i];
		}

		return temp;
	}
};

// Cont -> Container class
class Cont {

public:
	Matrix** mats;
	int matNums;

	Cont(int cellNums) {
		this->matNums = cellNums;
		mats = new Matrix*[cellNums];
	}

	void CreateCells() {
		for (int i = 0; i < matNums; i++)
		{
			mats[i] = new Matrix;
		}
	}

	~Cont(){

		for (int i = 0; i < matNums; i++)
			delete[] mats[i];

		delete[] mats;
	}
};