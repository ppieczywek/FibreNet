#pragma once
#include "Vector3.h"
/**
*	Structure holding  simulation settings defined in simulation input file.
**/
struct SimulationSettings
{
	int   types_num;
	int	  bc_type;
	double wall_contact_stiffness;
	double cutoff;
	double cutoff_sq;
	double box_unit_cell;
	double dt;
	double half_dt;
	double half_dt2;
	double lambda_dt;
	double gamma;
	double sigma;
	double damping;
	int x_cells_num;
	int y_cells_num;
	int z_cells_num;
	int xy_cells_num;
	double half_size_squared;
	Vector3 step;
	Vector3 size;
	Vector3 gravity;
};