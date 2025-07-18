/*This File Calculates the reciprocal energy (k!=0) using the integral method developed by Kawata et. al*/
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
const complex<double> t(0,1);

#define ENABLE_OMP 1

struct reciprocal_n_params{
  double* PosIons;
  double* ion_charges;
  int natoms;
  double betaa;
  double** box;
  int *K;
  double* Length;
};

double integrand_reciprocal(double h, void *params){
    // this integral is over the variable h;
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double *PosIons = p->PosIons;
    double *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    double **box = p->box;
    int *K = p->K;
    double* Length = p->Length;

    double reciprocal_energy_i=0;

    double G[3]={2*M_PI/Length[0], 2*M_PI/Length[1], h};
    double bet = 4*betaa*betaa;
    #if defined ENABLE_OMP
        #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i) collapse(2)
    #endif

    for (int k = -K[0]; k < K[0]+1; k++){
        for (int l = -K[1]; l < K[1]+1; l++){
            if((k==0) && (l==0))continue;
            complex<double> s_kh=0;
            for (int  i = 0; i < natoms; i++){
                double G_dot_r=k*G[0]*PosIons[3*i]+l*G[1]*PosIons[3*i+1]+G[2]*PosIons[3*i+2];
                s_kh+=ion_charges[i]*(cos(G_dot_r)+t*sin(G_dot_r));
            }
            double norm_sg = norm(s_kh);
            double norm_g = k*k*G[0]*G[0]+l*l*G[1]*G[1]+G[2]*G[2];
            reciprocal_energy_i+=norm_sg/(norm_g*exp(norm_g/(bet)));
        }
    }
    return reciprocal_energy_i;
}

double reciprocal_kawata(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* K) {
    // this is for Ui
    //Length of the sides of the unit cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10);

    gsl_function F;
    F.function = &integrand_reciprocal; // Set the function to integrate
    reciprocal_n_params params = {PosIons, ion_charges, natoms, betaa, box, K, Length};
    F.params = &params;

    double result, error;
    gsl_integration_qagi(&F, 1e-4, 1e-2, 10, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory
    result*=1/(Length[0]*Length[1]);
    return result;
}