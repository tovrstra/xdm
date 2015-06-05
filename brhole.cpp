

#include <cmath>
#include <stdexcept>


// 2/3*pi**(2/3)
#define CONST1 1.4300195980740170


inline double helper(double x, double rhs) {
    return x*exp(-2.0/3.0*x) - rhs*(x - 2);
}


void solve_brhole_isotropic(double* mgga, double* brhole, double threshold) {
    double rho = mgga[0];

    if (rho < threshold) {
        brhole[0] = 0.0;
        brhole[1] = 0.0;
        brhole[2] = 0.0;
        return;
    }

    double gradsq = mgga[1]*mgga[1] + mgga[2]*mgga[2] + mgga[3]*mgga[3];
    double d = 2.0*mgga[5] - gradsq/rho/4.0;
    double q = (mgga[4] - 2*d)/6.0;
    double rhs = CONST1*pow(rho, 5.0/3.0)/q;

    double x0 = 0.0;
    double f0 = helper(x0, rhs);
    double x1 = 1.0;
    double f1 = helper(x1, rhs);

    // Set the initial bracket.
    while (! ((f0 > 0) ^ (f1 > 0))) {
        x1 *= 2;
        f1 = helper(x1, rhs);
        if (x1 > 1e10) {
            throw std::overflow_error("Failed to find reasonable initial bracket.");
        }
    }

    // Make sure f0 is positive.
    if (f0 < 0) {
        double tmp = x0;
        x0 = x1;
        x1 = tmp;
        tmp = f0;
        f0 = f1;
        f1 = tmp;
    }

    // Refine the bracket until it is tight enough.
    int irep = 0;
    while (fabs(x0 - x1) > 1e-14) {
        // Select a new point in the bracket to become the new edge.
        double af0 = fabs(f0);
        double af1 = fabs(f1);
        double x2 = 0.0;
        if ((af0*10 < af1) | (af1*10 < af0)) {
            // Make the linear interpolation more robust and make a dumb esimate
            // of the root.
            x2 = (x0+x1)/2;
        } else {
            // Use linear interpolation to estimate the root.
            x2 = x0 - f0*(x1 - x0)/(f1 - f0);
        }
        // Evaluate f2 and decide which edge of the bracked to replace.
        double f2 = helper(x2, rhs);
        if (f2 > 0) {
            x0 = x2;
            f0 = f2;
        } else {
            x1 = x2;
            f1 = f2;
        }
        // throw an exception when there are more than 1000 iterations;
        if (irep > 1000)
            throw std::runtime_error("Failed to refine bracket in 1000 iterations.");
        irep++;
    }

    double x = (x0 + x1)/2;
    double b = 0.5*x*exp(-x/3.0)*pow(M_PI*rho, -1.0/3.0);
    double a = x/b;
    double A = pow(a, 3.0)/(8*M_PI);

    brhole[0] = A;
    brhole[1] = a;
    brhole[2] = b;
}
