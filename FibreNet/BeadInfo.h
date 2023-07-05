#pragma once
#include "Vector3.h"
/**
*	A strcuture holding additional data for bead structure.
**/
struct BeadInfo
{
	int			type;
	int			tag;
	int			loc_id;
	double		mass;
	double		mass_inv;
	double		lambda_dt_mass_inv;
	double		half_dt_mass_inv;
	double		half_dt2_mass_inv;
	double		radius;
	double		radius_sq;
	double		skin_radius;
	double		skin_radius_sq;
};