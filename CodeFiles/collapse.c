#pragma omp parallel for reduction(+: reci_energy) collapse(3)
for (int i = -K_x; i < K_x+1; i++){
    for (int j = -K_y; j< K_y+1; j++){
        for (int k = -K_z; k < K_z+1; k++){
            reci_energy+= ....;  // update energy
        }
    }
}