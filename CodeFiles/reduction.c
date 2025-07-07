#pragma omp parallel for reduction(+: real_energy)
for (int i = 1; i < natoms; i++){
    for (int j = 0; j < i; j++){
        ...
        ...
        real_energy+= ...; // add real space energy interaction between ith and jth atom
    }
}