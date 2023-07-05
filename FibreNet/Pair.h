#pragma once
/**
*	Structure holding  paramters of non-bondedn interactions.
**/
struct Pair
{
	int type;	// type of interactions model
	int p1;		// type of interacting bead
	int p2;		// type of interacting bead
	double c1;	// model coefficient
	double c2;  // model coefficient
	double c3;  // model coefficient
	double c4;  // model coefficient
	double c5;  // model coefficient
	double c6;  // model coefficient
};
