/*This file calculates the reciprocal space energy using the new method*/
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

#define REAL 0
#define IMAG 1
// Disable this declaration if openmp parallelization is not required, would not be helpful for smaller systems
#define ENABLE_OMP 1

double PM2DEwald(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* Grid, int *M, int* n){
    // n: order of b-spline interpolation
    // initializing the new variables for charge interpolation
    double *u, *x_direc, *y_direc, *z_direc;
    u= new double [natoms*3];
    x_direc= new double [natoms*Grid[0]];
    y_direc= new double [natoms*Grid[1]];
    z_direc= new double [natoms*Grid[2]];

    double L1 = box[0][0];
    double L2 = box[1][1];
    double L3 = box[2][2];
    int n_max=1;

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable using standard fftw syntax

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid[0]*Grid[1]*Grid[2]);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid[0]*Grid[1]*Grid[2]);
    fftw_plan p;
    p = fftw_plan_dft_3d(Grid[0],Grid[1],Grid[2], in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    #if defined ENABLE_OMP
        #pragma omp parallel for
    #endif
    // Calculating the fractional coordinates
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 3; j++){
            double x[3] = {PosIons[3*i],PosIons[3*i+1],PosIons[3*i+2]};
            u[3*i+j]=Grid[j]*dotProduct(x,G[j],3);
        }
    }

    #if defined ENABLE_OMP
        #pragma omp parallel for
    #endif
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < Grid[0]; k1++){
            x_direc[Grid[0]*i+k1]=0;
            for (int  n1 = -n_max; n1 < n_max+1; n1++){
                x_direc[Grid[0]*i+k1]+=M_n(u[3*i+0]-k1-n1*Grid[0],n[0]);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < Grid[1]; k2++){
            y_direc[Grid[1]*i+k2]=0;
            for (int  n2 = -n_max; n2 < n_max+1; n2++){
                y_direc[Grid[1]*i+k2]+=M_n(u[3*i+1]-k2-n2*Grid[1],n[1]);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < Grid[2]; k3++){
            z_direc[Grid[2]*i+k3]=0;
            for (int  n3 = -n_max; n3 < n_max+1; n3++){
                z_direc[Grid[2]*i+k3]+=M_n(u[3*i+2]-k3-n3*Grid[2],n[2]);
            }
        }
    }

    // initializing the "in" vector with zero values
    for (int tx = 0; tx < Grid[0]; tx++){
        for (int ty = 0; ty < Grid[1]; ty++){
            for (int tz = 0; tz < Grid[2]; tz++){
                in[tx * (Grid[2] * Grid[1]) + ty * Grid[2] + tz][0] = 0.0;
            }
        }
    }

    #if defined ENABLE_OMP
        #pragma omp parallel for 
    #endif
    // Final Q Matrix
    for (int j = 0; j < natoms; j++){
        if (ion_charges[j] == 0)continue;
        for (int tx = 0; tx < Grid[0]; tx++){
            if (x_direc[Grid[0]*j+tx] == 0)continue;

            for (int ty = 0; ty < Grid[1]; ty++){
                if (y_direc[Grid[1]*j+ty] == 0)continue;

                for (int tz = 0; tz < Grid[2]; tz++){
                    if (z_direc[Grid[2]*j+tz] == 0)continue;
                    #if defined ENABLE_OMP
                        #pragma omp atomic update
                    #endif
                    in[tx * (Grid[2] * Grid[1]) + ty * Grid[2] + tz][REAL] += ion_charges[j] * x_direc[Grid[0]*j+tx] * y_direc[Grid[1]*j+ty] * z_direc[Grid[2]*j+tz];
                }
            }
        }
    }

    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup();

    double energy=0;
    // collapse doesn't makes a difference here much; dynamic and runtime give the same time 
    int ii,jj,kk;
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: energy) collapse(3)
    #endif
    for (int i = -M[0]; i < M[0]+1; i++){
        for (int j = -M[1]; j< M[1]+1; j++){
            for (int k = -M[2]; k < M[2]+1; k++){
                if(i==0&&j==0&&k==0)continue;
                int ic,jc,kc;
                if(i<0) {ii=Grid[0]+i;ic=(2*M[0]+1)+i;}
                else {ii=i;ic=i;}
                if(j<0) {jj=Grid[1]+j;jc=(2*M[1]+1)+j;}
                else  {jj=j;jc=j;}
                if(k<0) {kk=Grid[2]+k;kc=(2*M[2]+1)+k;}
                else  {kk=k;kc=k;}

                int temp=ii * (Grid[2] * Grid[1]) + jj * Grid[2] + kk;
                int tempexpfactor = ic * ((2*M[1]+1) * (2*M[2]+1)) + jc * (2*M[2]+1) + kc;

                double norm_FQ=out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];

                energy += norm_FQ*ExpFactorInterpolated[tempexpfactor]; //update energy
            }
        }
    }
    return energy*sqrt(M_PI)/(L1*L2);
}