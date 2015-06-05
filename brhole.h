

#ifndef BRHOLE_H
#define BRHOLE_H

/**
    @brief
        Solver for the isotropic Becke-Rousel hole model.

    @param mgga
        A pointer to six doubles with the MGGA inpuyt data for one grid point:
        rho, drho_dx, drho_dy, drho_dz, lapl, tau.

    @param brhole
        A pointer to three double that will contain the output arguments for
        the Becke-Rousel hole model: A, a, b.

    @param threshold
        If the density falls below this threshold, all Becke-Rousel parameters
        are set to zero and no attempt is made to solve the equations.
*/
void solve_brhole_isotropic(double* mgga, double* brhole, double threshold);

#endif
