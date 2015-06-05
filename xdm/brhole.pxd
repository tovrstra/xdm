

cdef extern from "brhole.h":
    void solve_brhole_isotropic(double* mgga, double* brhole, double threshold)
