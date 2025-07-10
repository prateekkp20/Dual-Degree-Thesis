#include<omp.h>
#include<stdio.h>

int main() {
    omp_set_num_threads(4);  // Set number of threads to 4
    const int N = 100;
    double result[N];

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        result[i] = i * 2.0;  
    }
    return 0;
}