/*This Method is for the correction methods that we have introduced using the tophat function for ewald summation without the use of any particle-mesh methods*/
double reciprocal_modified(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int *K){
    double reci_energy=0;
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: reci_energy) collapse(3)
    #endif
    for (int kx = -K[0]; kx < K[0]+1; kx++){
        for (int ky = -K[1]; ky < K[1]+1; ky++){
            for (int kz = -K[2]; kz < K[2]+1; kz++){
                if((kx==0) && (ky==0) && (kz==0))continue;

                int ii,jj,kk;
                if(kx<0) ii=(2*K[0]+1)+kx;
                else ii=kx;
                if(ky<0) jj=(2*K[1]+1)+ky;
                else  jj=ky;
                if(kz<0) kk=(2*K[2]+1)+kz;
                else  kk=kz;
                int temp=ii * ((2*K[2]+1) * (2*K[1]+1)) + jj * (2*K[2]+1) + kk;

                complex<double> sg=0;
                for (int  i = 0; i < natoms; i++){
                    double G_dot_r=2*M_PI*(kx*G[0][0]*PosIons[3*i]+ky*G[1][1]*PosIons[3*i+1]+kz*G[2][2]*PosIons[3*i+2]);
                    sg+=ion_charges[i]*(cos(G_dot_r)+t*sin(G_dot_r));
                }
                double norm_sg = norm(sg);
                
                //update energy
                reci_energy+=ExpFactor[temp]*norm_sg;
            }
        }
    }
    double L1 = box[0][0];
    double L2 = box[1][1];
    double L3 = box[2][2];
    return reci_energy*sqrt(M_PI)/(L1*L2);
}