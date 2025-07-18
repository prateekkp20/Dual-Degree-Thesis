//***********************************************
// function to calculate distances according to minimum image convention for a 2D symmetric system
//***********************************************

#include "libinclude.h" 

double dist(double* __restrict__ PosIons2, int atom1, int atom2, double** __restrict__ box){

	double Dx, Dy, Dz;
	double Dx1, Dy1, Dz1;
    double distatoms;

	Dx = PosIons2[atom1*3] - PosIons2[atom2*3];
	Dy = PosIons2[atom1*3+1] - PosIons2[atom2*3+1];
	Dz = PosIons2[atom1*3+2] - PosIons2[atom2*3+2];

	double a  = ceil(Dx/box[0][0]-0.5), b = ceil(Dy/box[1][1]-0.5);

    Dx1 = Dx - box[0][0]*a - box[1][0]*b;
	Dy1 = Dy - box[0][1]*a - box[1][1]*b;
	Dz1 = Dz - box[0][2]*a - box[1][2]*b;

	distatoms = sqrt(Dx1*Dx1 + Dy1*Dy1 + Dz1*Dz1);
	return distatoms;
	
}