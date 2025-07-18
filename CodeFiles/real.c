function realExact(Positions, Charges, Natoms, alpha, box, cutoff):
    
    for i from 1 to Natoms:
        for j from 1 to i:
            R = dist(Positions, i, j ,box)
            
            // skip if R > cutoff

            real_energy + = charge_i*charge_j*erfc(alpha*R)/R
    
    return real_energy