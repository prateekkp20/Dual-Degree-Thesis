// this code implements the reciprocal (k!=0) energy using the bspline method with 2D Fourier Transform and 1D Fourier Integral
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

#define ENABLE_OMP 1
const complex<double> t(0.0,1.0);
complex<double> Coeff_H_nz;
double FourPiPi = 4*M_PI*M_PI;
double deno;
double *ReciVector;

struct reciprocal_n_params {
    double* ion_charges;
    int natoms;
    int *K;
    int *Grid;
    int *n;
    double **x_direc, **y_direc, **z_direc;
    double * TZ ;
    int GridZ;
};

double reciprocal_ft_integrand(double h, void *params){
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    int *K = p->K;
    int *Grid = p->Grid;
    // n: order of b-spline interpolation
    int *n = p->n;
    auto x_direc=p->x_direc,y_direc=p->y_direc,z_direc=p->z_direc;
    double *TZ=p->TZ;
    int GridZ = p->GridZ;

    double reciprocal_energy_i=0;

    // the fourier integral of z_direc vector for every ith atom
    complex<double>* fz_i_h = new complex<double> [natoms];
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int  i = 0; i < natoms; i++){
        for (int  tz = 0; tz < GridZ; tz++){
            fz_i_h[i]+=z_direc[i][tz]*exp(h*TZ[tz]*t);
        }
    }

    fftw_complex *in;
    fftw_complex *out;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid[0]*Grid[1]);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid[0]*Grid[1]);

    fftw_plan plan;
    plan = fftw_plan_dft_2d(Grid[0], Grid[1], in ,out, FFTW_BACKWARD, FFTW_ESTIMATE);

    complex<double> temp;
    #if defined ENABLE_OMP
        #pragma omp parallel for private(temp)
    #endif
    for (int i = 0; i < natoms; i++){
        for (int tx = 0; tx < Grid[0]; tx++){
            if(x_direc[i][tx]==0)continue;
            for (int ty = 0; ty < Grid[1]; ty++){   
                if(y_direc[i][ty]==0)continue;
                temp = ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i];
                #if defined ENABLE_OMP
                    #pragma omp atomic update
                #endif
                    in[Grid[1]*tx+ty][0]+=temp.real();
                #if defined ENABLE_OMP
                    #pragma omp atomic update
                #endif
                    in[Grid[1]*tx+ty][1]+=temp.imag();
            }
        }
    }

    fftw_execute(plan);

    Coeff_H_nz = Coeff(h,8);
    #pragma omp parallel for schedule(runtime) reduction(+: reciprocal_energy_i) collapse(2)
    for (int i = -K[0]; i < K[0]+1; i++){
        for (int j = -K[1]; j< K[1]+1; j++){
            int ii,ic,jj,jc;
            if(i==0&&j==0)continue;
            if(i<0){ii=Grid[0]+i;ic=(2*K[0]+1)+i;}
            else {ii=i;ic=i;}
            if(j<0){jj=Grid[1]+j;jc=(2*K[1]+1)+j;}
            else {jj=j;jc=j;}
            int temp = Grid[1]*ii+jj;
            double factor = ReciVector[ic*(2*K[1]+1)+jc] + h*h;
            double norm_FQ = norm(out[temp][0] + t*out[temp][1]);
            reciprocal_energy_i+= norm_FQ  * norm(CoeffX[ic]*CoeffY[jc]*Coeff_H_nz) / (factor*exp(factor/deno));
        }
    }
    return reciprocal_energy_i;
}

// Main Function to Calculate the Reciprocal Energy (k!=0)
double reciprocal_fft(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* K, int *Grid, int *n){

    // Edge lengths of the cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};

    deno = 4*betaa*betaa;

    int ii,ic,jj,jc;
    ReciVector = new double [(2*K[0]+1)*(2*K[1]+1)];
    for (int i = -K[0]; i < K[0]+1; i++){
        for (int j = -K[1]; j< K[1]+1; j++){
                if(i<0){ic=(2*K[0]+1)+i;}
                else {ic=i;}
                if(j<0){jc=(2*K[1]+1)+j;}
                else {jc=j;}
                ReciVector[ic*(2*K[1]+1)+jc] = FourPiPi * (i*i*G[0][0]*G[0][0]+j*j*G[1][1]*G[1][1]);
        }
    }

    // initializing the new variables
    // u: the fractional coordinates in x and y directions
    // x_direc, y_direc, z_direc: the cofficients in the x,y and z directions for the Q Matrix
    int GridZ = Length[2]+n[2]+1;
    double **u,**x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];
    z_direc = new double * [natoms]; 
    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2]; // We only need these in x and y direction 
        x_direc[i] = new double  [Grid[0]];
        y_direc[i] = new double  [Grid[1]];
        z_direc[i] = new double  [GridZ]; // tz varies from -n to Zmax(Lz)
    }

    // Calculating the fractional coordinates in x and y directions
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            double x[3] = {PosIons[3*i],PosIons[3*i+1],PosIons[3*i+2]};
            u[i][j]=Grid[j]*dotProduct(x,G[j],3);
        }
    }
    double * TZ = linspace(-n[2],(int)Length[2],1);
    int l_max=1;

    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  tx = 0; tx < Grid[0]; tx++){
            x_direc[i][tx]=0;
            for (int  lx = -l_max; lx < l_max+1; lx++){
                x_direc[i][tx]+=M_n(u[i][0]-tx-lx*Grid[0],n[0]);
            }
        }
        // for Y direction
        for (int  ty = 0; ty < Grid[1]; ty++){
            y_direc[i][ty]=0;
            for (int  ly = -l_max; ly < l_max+1; ly++){
                y_direc[i][ty]+=M_n(u[i][1]-ty-ly*Grid[1],n[1]);
            }
        }
        // for Z direction
        for (int  tz = 0; tz < GridZ; tz++){
            z_direc[i][tz]=M_n(PosIons[3*i+2]-TZ[tz],n[2]);
        }
    }
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &reciprocal_ft_integrand; // Set the function to integrate
    reciprocal_n_params params = {ion_charges, natoms, K, Grid, n, x_direc, y_direc, z_direc, TZ, GridZ};
    F.params = &params;
    double result, error;
    gsl_integration_qagi(&F, 1e-4, 1e-2, 200, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory
    
    delete [] TZ;
    for (int  i = 0; i < natoms; i++){
        delete [] u[i];
        delete [] x_direc[i];
        delete [] y_direc[i];
        delete [] z_direc[i];
    }
    delete [] u;
    delete [] x_direc;
    delete [] y_direc;
    delete [] z_direc;
    delete [] ReciVector;

    result*=1/(Length[0]*Length[1]);
    return result;
}