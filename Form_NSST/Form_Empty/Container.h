#pragma once
#include <stdlib.h>

class Matrix {

public:
	int height;
	int width;
	int depth;

	double* mat;

	void CreateMatrix(int height, int width, int depth = 1) {

		this->height = height;
		this->width = width;
		this->depth = depth;
		mat = (double*)calloc(height * width * depth, sizeof(double));
	}

	int GetSize() { return height * width; }

	double operator()(int row, int col) {
		return this->mat[row * width + col];
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