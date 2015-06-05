#cython: embedsignature=True
'''C++ extensions'''

import numpy as np
cimport numpy as np
np.import_array()

cimport brhole


__all__ = ['solve_brhole_isotropic', 'solve_brhole_isotropic_array']


def solve_brhole_isotropic(np.ndarray[double, ndim=1] mgga not None,
                           double threshold=1e-6):
    assert mgga.shape[0] == 6
    assert mgga.flags['C_CONTIGUOUS']
    cdef np.ndarray[double, ndim=1] output = np.zeros(3, float)
    brhole.solve_brhole_isotropic(&mgga[0], &output[0], threshold)
    return output


def solve_brhole_isotropic_array(np.ndarray[double, ndim=2] mgga not None,
                           double threshold=1e-6):
    cdef long npoint = mgga.shape[0]
    assert mgga.shape[1] == 6
    assert mgga.flags['C_CONTIGUOUS']
    cdef np.ndarray[double, ndim=2] output = np.zeros((npoint, 3), float)

    cdef long ipoint = 0
    for ipoint in xrange(npoint):
        brhole.solve_brhole_isotropic(&mgga[ipoint,0], &output[ipoint,0], threshold)
    return output
