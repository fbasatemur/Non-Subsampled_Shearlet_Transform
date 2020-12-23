#pragma once

class Cell {

public:
	double* matx;
	int cols;
	int rows;
	int depth;
	int matxSize;
	int matxCount;

	Cell();
	Cell(int width, int height, int depth = 1);
	~Cell();

	double* CreateMatrix(int height, int width, int depth = 1);
	int getSize();

	Cell* operator*(double scalarValue);

};

//Non-Member Function
//Allocate the new cell pointer with this function.
Cell* newCell(int matxCount)
{
	Cell* cell = new Cell[matxCount];
	cell->matxCount = matxCount;
	return cell;
}