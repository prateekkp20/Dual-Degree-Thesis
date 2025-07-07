#pragma omp parallel for
for (int j = 0; j < natoms; ++j) {
    for (int tx = 0; tx < Grid[0]; ++tx) {
        for (int ty = 0; ty < Grid[1]; ++ty) {
            for (int tz = 0; tz < Grid[2]; ++tz) {
                ...
                double w = compute_weight(j, tx, ty, tz);
                // Safely accumulate weight onto the grid point using atomic
                #pragma omp atomic update
                Q[tx][ty][tz] += w;
            }
        }
    }
}