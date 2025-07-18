function dist(Positions, atom1, atom2, box):
    // Compute distance between atom1 and atom2
    Dx = ...
    Dy = ...
    Dz = ...

    // Apply Minimum Image Convention (MIC)
    // General MIC formula for 1D:
    //     delta_x = delta_x - L * round(delta_x / L);

    
    Dx1 = Dx - L1 - ... // ceil() function is used to round
    Dy1 = Dy - L2 - ...
    Dz1 = Dz - L3 - ...

    // Final Euclidean distance after applying MIC
    distance = sqrt(Dx1^2 + Dy1^2 + Dz1^2)
    return distance
