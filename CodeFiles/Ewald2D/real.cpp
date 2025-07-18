#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double real(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
    double real_energy=0;

    #if defined ENABLE_OMP
        #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
    #endif

        for (int i = 1; i < natoms; i++){
            for (int j = 0; j < i; j++){
                    double modR=dist(PosIons,i,j,box);
                    if(modR>cutoff)continue;

                    // Erfc Approximations
                    double val = betaa*modR;
                    double exp_x2 = exp(-val*val);
                    double t, t1 =  t  = 1/(1+0.3275911*val);
                    double erfcx = exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
                    real_energy+=(ion_charges[i]*ion_charges[j]*erfcx)/modR;
            }
        }
    
    return real_energy;
}

double realExact(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
    double real_energy=0;

    #if defined ENABLE_OMP
        #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
    #endif

        for (int i = 1; i < natoms; i++){
            for (int j = 0; j < i; j++){
                    double modR=dist(PosIons,i,j,box);
                    if(modR>cutoff)continue;

                    real_energy+=(ion_charges[i]*ion_charges[j])*erfc(betaa*modR)/modR;
            }
        }
    
    return real_energy;
}