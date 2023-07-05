#pragma once
/**
*	A strcuture holding data of if single angular bond.
**/
struct AngleBond
{
	int						type;			//	bond type 
	int						tag;			//	
	int						status;			//	status of a bond when destructible type is used 
	int						b1;				//	first bead id
	int						b2;				//	second bead id
	int						b3;				//	third bead id
	int						s1;				//	first spring id
	int						s2;				//	second spring id
	double					angle;			//	bond rest angle
	double					current_angle;	//	bond currnet angle
	double					c1;				//	bond stiffness
	double					c2;				//
	double					c3;				//	
};