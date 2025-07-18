#include "libinclude.h"
#include "const.h"
#include "fundec.h"

// Function to calculate the dot product of two vectors of size n
template<typename T>
double dotProduct(T v1, T v2, size_t n){
    double result = 0;
    for (size_t i = 0; i < n; ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}
template double dotProduct<float*>(float* v1, float* v2, size_t size);
template double dotProduct<double*>(double* v1, double* v2, size_t size);

// For the reciprocal space energy (k=0)
double F_0(double val){
    double a = 1-exp(-val*val);
    double b = sqrt(M_PI)*val*erf(val);
    return a-b;
}

// Function to calculate the cross product of two vectors and store it in the "out" vector
template<typename T1>
void crossProduct(T1 v_A, T1 v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}
template void crossProduct<float*>(float*, float*, double*);
template void crossProduct<double*>(double*, double*, double*);

/**
 * @brief Computes the value of the B-spline basis function of degree (n-1) at position u using the Cox-de Boor recursion formula.
 *
 * This function recursively evaluates the B-spline basis function M_n(u, n) for a given position u and order n.
 * The implementation follows the Cox-de Boor recursion formula, which is fundamental in the construction of B-spline curves.
 *
 * @param u The position at which to evaluate the B-spline basis function.
 * @param n The order of the B-spline (degree = n-1). Must be >= 2.
 * @return The value of the B-spline basis function at position u. Returns 0 if u is outside the valid range [0, n].
 *
 * @note For n == 2, the function returns the linear B-spline (tent function).
 * @note For n < 2, the function returns 0 as the basis is not defined.
 */
double M_n(double u, int n){
    if(n<2)return 0;
    else if(n==2){
        if(u<0 || u>n) return 0;
        else{
            return 1-abs(u-1);
        }
    }
    else{
        if(u<0 || u>n) return 0;
        else{
            return (u*M_n(u,n-1)/(n-1))+((n-u)*M_n(u-1,n-1)/(n-1));
        }
    }
}


/**
 * @brief Computes B-spline coefficients for Ewald summation.
 *
 * This function calculates the B-spline coefficient for given parameters v and w,
 *
 * @param v The phase or frequency parameter.
 * @param w The order or degree of the B-spline.
 * @return The computed B-spline coefficient as a complex number.
 */
complex<double>Coeff(double v, double w){
    const complex<double> t(0, 1);
    complex<double> bi_mi=exp(t*v*(w-1));
    complex<double> denox;
    for (double f = 0; f < w-1; f++){
        denox+=M_n(f+1,w)*exp(t*f*v);
    }
    if (denox == complex<double>(0, 0)) return 0;
    bi_mi/=denox;
    return bi_mi;
}