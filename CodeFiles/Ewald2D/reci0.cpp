/*This file calculates the reciprocal space energy (k=0), method given by parry*/
#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double reci0(double *PosIons2, double *ion_charges, int natoms, double betaa, double **box){
    double energy = 0;
    double Length[3]={box[0][0],box[1][1],box[2][2]};

    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: energy)
    #endif
    
    for (int  i = 1; i < natoms; i++){
        for (int j = 0; j < i; j++){
            energy+=ion_charges[i]*ion_charges[j]*F_0((PosIons2[i*3+2]-PosIons2[j*3+2])*betaa);
        }
    }
    
    return 2*sqrt(M_PI)*energy/(betaa*Length[0]*Length[1]);
}