#include "Cell.h"

Cell::Cell()
{
	//ctor
}

Cell::Cell(int width, int height, int depth) {
	this->cols = width;
	this->rows = height;
	this->depth = depth;
	this->matxSize = cols * rows * depth;
	this->matx = new double[matxSize];
}

Cell::~Cell()
{
	delete[] matx;
}

double* Cell::CreateMatrix(int height, int width, int depth) {
	this->cols = width;
	this->rows = height;
	this->depth = depth;
	this->matxSize = cols * rows * depth;
	this->matx = new double[matxSize];
	return this->matx;
}

int Cell::getSize() {
	return this->cols * this->rows;
}

Cell* Cell::operator* (double value)
{
	for (int i = 0; i < this->matxSize; i++)
		this->matx[i] *= value;
	return this;
}
