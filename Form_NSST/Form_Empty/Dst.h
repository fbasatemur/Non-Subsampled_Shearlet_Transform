#pragma once

class Dst {

public:
	double** mats;
	int* matsDepth;
	int cols;
	int rows;
	int matCount = 0;
	
	Dst(int width, int height, int matrixCount) {
		this->cols = width;
		this->rows = height;
		this->matCount = matrixCount;

		this->mats = new double* [matrixCount];
		this->matsDepth = new int[matrixCount];
	}

	~Dst(){
		delete[] matsDepth;

		for (int i = 0; i < matCount; i++)
			delete[] mats[i];
	}

	void setMatDepth(int matIndex, int depth) {
		this->matsDepth[matIndex] = depth;
	}

	int getMatDepth(int matIndex) {
		return this->matsDepth[matIndex];
	}

	int getSize() {
		return this->cols* this->rows;
	}

};
