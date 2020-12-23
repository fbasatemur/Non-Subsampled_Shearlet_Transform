#pragma once

struct Matrix {

	int height;
	int width;
	int depth;

	double* mat;

	void CreateMatrix(int height, int width, int depth = 1) {

		this->height = height;
		this->width = width;
		this->depth = depth;
		mat = new double[height * width * depth];
	}

	int GetSize() { return height * width; }
};

// Cont -> Container class
class Cont {

public:
	Matrix* mats;
	int matNums;

	Cont(int cellNums) {
		this->matNums = cellNums;
		mats = new Matrix[cellNums];
	}

	~Cont(){

		for (int i = 0; i < matNums; i++)
			delete[] mats[i].mat;

		delete[] mats;
	}

};