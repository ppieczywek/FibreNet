#pragma once
#include "Vector3.h"
/**
*	A structure holding addresses of surrounding cells for each of collision grid cells.
**/
struct GridCellNh
{
	int n[13];
};

/**
*	Single bead collision bucket, holding ID's of petentially colliding beads.
**/
struct ParticleNh
{
	int n[256];
};

/**
*	Single collision grid cell bucket, holding ID's of beads occupying this cell.
**/
struct GridCell
{
	int n[64];
};