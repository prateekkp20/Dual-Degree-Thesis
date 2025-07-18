// This file has the functions, and data structures to calculate the Screening Factor for our modified method 
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

struct C_params{
    double t; double L; double gamma; double n;
};

struct constant_params{
    int kx; int ky; int kz; double lx; double ly; double lz; double gamma;
};

double tophat(double s, double gamma, double L){
    double result = 1/(1+exp(-gamma*(0.5*L+s))) + 1/(1+exp(-gamma*(0.5*L-s))) - 1; // sigmoid type function
    // double result = tanh(gamma*(s+0.5*L)) - tanh(gamma*(s-0.5*L));
    return result;
}

double C_integrand(double s, void *params){
    C_params* p = (struct C_params*)params;
    double t = p->t;
    double L = p->L;
    double gamma = p->gamma;
    double n = p->n;

    double result=exp(-s*t*s*t)*cos(2*M_PI*n*s/L)*tophat(s,gamma,L);
    return result;
}

double C(double t, double n, double gamma, double L){

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &C_integrand; // Set the function to integrate
    C_params params = {t, L, gamma, n};
    F.params = &params;
    double result, error;
    gsl_integration_qagi(&F, 1e-8, 1e-6, 200, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result/L;
}

double constant_integrand(double t, void *params){
    constant_params* p = (struct constant_params*)params;
    int kx = p->kx; int ky = p->ky; int kz = p->kz; double lx = p->lx; double ly = p->ly; double lz = p->lz; double gamma = p->gamma;
    double m[3];
    for (int t = 0; t < 3; t++){
        m[t]=kx*G[0][t]+ky*G[1][t];   
    }
    double g=M_PI*M_PI*dotProduct(m,m,3);
    double result = (1/(t*t))*C(t,kz,gamma,lz)*exp(-g/(t*t));
    return result;
}

double constantterm(int kx, int ky, int kz, double lx, double ly, double lz, double beta, double gamma){

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &constant_integrand; // Set the function to integrate
    constant_params params = {kx, ky, kz, lx, ly, lz, gamma};
    F.params = &params;
    double result, error;
    size_t neval;
    gsl_integration_qng(&F, 0, beta, 1e-4, 1e-4, &result, &error, &neval);
    // int gsl_integration_qng(const gsl_function *f, double a, double b, double epsabs, double epsrel, double *result, double *abserr, size_t *neval)

    // gsl_integration_qag(&F,0, beta, 1e-4, 1e-2, 200, 3, workspace, &result, &error);// this has less error 
    // int gsl_integration_qag(const gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr)

    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result;
}

/*
 * @brief Determines the minimum vacuum size required to fully contain the simulation box within the tophat window.
 *
 * This function uses a binary search algorithm to find the smallest value of the tophat window length (vacuum)
 * such that the simulation box of given side length lies entirely inside the window, as defined by the tophat function.
 * The function checks if the initial maximum vacuum is sufficient to contain the box; if not, it exits with an error.
 * The search is performed between twice the side length and the provided maximum vacuum, with the result rounded to
 * the nearest 0.1 and an optional margin added for safety.
 *
 * @param SideLength The length of the simulation box side to be contained within the tophat window.
 * @param gamma      The sharpness parameter for the tophat (sigmoid) function.
 * @param maxVacuum  The maximum allowed vacuum size to consider during the search.
 * @param margin     An additional margin to add to the computed vacuum size for safety.
 * @return           The minimum vacuum size (plus margin) required to fully contain the simulation box.
 *
 * @note The function assumes that the tophat function returns 1 when the box is fully contained.
 *       If the initial maxVacuum is insufficient, the function prints an error and exits.
 */
double vacuum(double SideLength, double gamma, double maxVacuum, double margin){
    if(tophat(SideLength,gamma,maxVacuum)!=1){
        printf("Error: Simulation box does not fit inside the tophat window. Exiting.\n");
        exit(EXIT_FAILURE);
    }
    double high = maxVacuum, low = 2*SideLength;
    while(high - low > 1){
        double mid = (high+low)/2;
        if(tophat(SideLength,gamma,mid) == 1){
            high = mid;
        }
        else{
            low = mid;
        }
    }
    return margin + round(high * 10) / 10.0;
}

/*Function for configuration independent computations for the reciprocal space summation*/
/* Exp(-|G|)/|G| or the screening function term in the reciprocal loop for direct ewald*/
// The negative indices are stored from the end of the array
void ScreenFunction(int *Kvec, double gamma, double beta, double** boxcell){
    
    // side lengths
    double L1 = boxcell[0][0];
    double L2 = boxcell[1][1];
    double L3 = boxcell[2][2];
    
    #pragma omp parallel for schedule(runtime) collapse(3)
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        for (int j = -Kvec[1]; j< Kvec[1]+1; j++){
            for (int k = -Kvec[2]; k < Kvec[2]+1; k++){
                if(i==0&&j==0&&k==0)continue;
				int ii,jj,kk;
				if(i<0) ii=(2*Kvec[0]+1)+i;
                else ii=i;
                if(j<0) jj=(2*Kvec[1]+1)+j;
                else  jj=j;
                if(k<0) kk=(2*Kvec[2]+1)+k;
                else  kk=k;
				int temp=ii * ((2*Kvec[2]+1) * (2*Kvec[1]+1)) + jj * (2*Kvec[2]+1) + kk;
				ExpFactor[temp] = constantterm(i,j,k,L1,L2,L3,beta,gamma);
			}
		}
	}
    return;
}

/*
 @brief Computes and initializes the SPME (Smooth Particle Mesh Ewald) coefficients for each spatial dimension. Cofficients Bi[mi] of the bspline interpolation for euler exponential splines. Refer to Essmann et al.
 *
 * This function allocates and fills the coefficient arrays (CoeffX, CoeffY, CoeffZ) used in the SPME method
 * for Ewald summation in periodic systems. For each spatial direction (x, y, z), it computes the coefficients
 * based on the provided grid size, interpolation order, and the number of reciprocal lattice vectors (Kvec).
 * The coefficients are calculated using the Coeff function for each relevant wave vector.
 *
 * @param Kvec  Pointer to an array of 3 integers specifying the number of reciprocal lattice vectors
 *              in the x, y, and z directions, respectively.
 * @param grid  Pointer to an array of 3 integers specifying the grid size in the x, y, and z directions.
 * @param order Pointer to an array of 3 integers specifying the interpolation order for each direction.
 */
void SPME_Coeff(int *Kvec, int* grid, int* order){

	double TwoPi_Gridx = 2*M_PI/grid[0];
    double TwoPi_Gridy = 2*M_PI/grid[1];
    double TwoPi_Gridz = 2*M_PI/grid[2];
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[0]+1)+i;}
        else {ic=i;}
        CoeffX[ic] = Coeff(TwoPi_Gridx*i,order[0]);
    }
	for (int i = -Kvec[1]; i < Kvec[1]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[1]+1)+i;}
        else {ic=i;}
        CoeffY[ic] = Coeff(TwoPi_Gridy*i,order[1]);
    }
	for (int i = -Kvec[2]; i < Kvec[2]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[2]+1)+i;}
        else {ic=i;}
        CoeffZ[ic] = Coeff(TwoPi_Gridz*i,order[2]);
    }
}

/*Function for configuration independent computations for the reciprocal space summation*/
/* B(m1,m2,m3)*Exp(-|G|)/|G| term in the reciprocal loop for SPME*/
// The negative indices are stored from the end of the array
void ScreenFunctionSPME(int *Kvec, int* grid, int* order, double gamma, double beta, double** boxcell){

    // side lengths
    double L1 = boxcell[0][0];
    double L2 = boxcell[1][1];
    double L3 = boxcell[2][2];

    #pragma omp parallel for schedule(runtime) collapse(3)
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        for (int j = -Kvec[1]; j< Kvec[1]+1; j++){
            for (int k = -Kvec[2]; k < Kvec[2]+1; k++){
				if(i==0&&j==0&&k==0)continue;
				int ii,jj,kk;
				if(i<0) ii=(2*Kvec[0]+1)+i;
                else ii=i;
                if(j<0) jj=(2*Kvec[1]+1)+j;
                else  jj=j;
                if(k<0) kk=(2*Kvec[2]+1)+k;
                else  kk=k;
				int temp=ii * ((2*Kvec[2]+1) * (2*Kvec[1]+1)) + jj * (2*Kvec[2]+1) + kk;
				ExpFactorInterpolated[temp] = norm(CoeffX[ii]*CoeffY[jj]*CoeffZ[kk])*constantterm(i,j,k,L1,L2,L3,beta,gamma);
			}
		}
	}
    return;
}