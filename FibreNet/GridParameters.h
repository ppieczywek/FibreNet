#pragma once
#include "Vector3.h"

struct GridParameters
{
	int x_cells_num;
	int y_cells_num;
	int z_cells_num;
	int xy_cells_num;
	double half_size_squared;
	Vector3 step;
	Vector3 size;
};