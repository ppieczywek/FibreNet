#pragma once
#include "Vector3.h"
/**
*	A strcuture holding data for harmonic bond structure.
**/
struct SpringStructure
{	
	int						type;
	int						tag;
	int						status;
	int						p1;
	int						p2;
	double					rest_length;
	double					stiffness;
	double					damping;
	double					c1;
	double					c2;
};