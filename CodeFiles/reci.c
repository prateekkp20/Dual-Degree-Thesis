function PM2DEwald(Positions, Charges, natoms, alpha, box, grid, Kx, Ky, Kz, order)
    
    // Allocate complex arrays `in[0..grid[0]-1][0..grid[1]-1][0..grid[2]-1]' and `out[0..grid[0]-1][0..grid[1]-1][0..grid[2]-1]' initialized to zero
    // Create FFTW plan `p' for 3D fwd transform of `in', output to `out'

    // Store Fractional Coordinates in U array
    for i from 1 to natoms:
        // Set U <- grid . position_i

    // Calculate Q-matrix coefficients in each direction
    for i from 1 to natoms:
        for k1 from 1 to grid[0]:
            Qx_i_k1 += // Assign the total B-spline function, defined in Eq.(2.11)
        for k2 from 1 to grid[1]:
            Qy_i_k2 += ...
        for k3 from 1 to grid[2]:
            Qz_i_k3 += ...

    // Build the final Q-matrix (charge assignment)
    for i from 1 to natoms:
        for (k1,k2,k3) from (1,1,1) to grid: // Short-hand written for nested loops 
            in[tx,ty,tz] = charge_i * Qx_i_k1 * Qy_i_k2 * Qz_i_k3

    // Run FFT
    // Execute FFTW plan `p' on `in', store result in `out'
    // Destroy FFTW plan `p' and cleanup

    // Compute Reciprocal Space Energy
     for i from -Kx to Kx:
        for j from -Ky to Ky:
            for k from -Kz to Kz:
                continue if i=0 and j=0 and k=0
                ...
                energy+=ExpFactorInterpolated[i, j, k] * norm(out[i, j, k]) 
    
    return energy * sqrt(pi) / (L1*L2)