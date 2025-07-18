#include "libinclude.h"
#include "const.h"
#include "fundec.h"

double self(int n_atomtype, int *natoms_type, double *chargs, double betaa){
    double self_energy=0;
    #pragma omp SIMD
    for (int i = 0; i < n_atomtype; i++){
        self_energy+=natoms_type[i]*pow(chargs[i],2);
    }
    return self_energy*(-betaa)/sqrt(M_PI);
}