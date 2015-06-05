#!/usr/bin/env python

import sys
from horton import *
import h5py as h5, numpy as np

log.set_level(log.silent)

pold_unit = angstrom**3
pold_map = {
    1: 0.6668,
    2: 0.2051,
    3: 24.33,
    4: 5.6,
    5: 3.03,
    6: 1.76,
    7: 1.1,
    8: 0.802,
    9: 0.557,
    10: 0.3956,
    11: 24.11,
    12: 10.6,
    13: 6.8,
    14: 5.38,
    15: 3.63,
    16: 2.9,
    17: 2.18,
    18: 1.6411,
    19: 43.4,
    20: 22.8,
    21: 17.8,
    22: 14.6,
    23: 12.4,
    24: 11.6,
    25: 9.4,
    26: 8.4,
    27: 7.5,
    28: 6.8,
    29: 6.2,
    30: 5.75,
    31: 8.12,
    32: 6.07,
    33: 4.31,
    34: 3.77,
    35: 3.05,
    36: 2.4844,
    37: 47.3,
    38: 27.6,
    39: 22.7,
    40: 17.9,
    41: 15.7,
    42: 12.8,
    43: 11.4,
    44: 9.6,
    45: 8.6,
    46: 4.8,
    47: 7.2,
    48: 7.36,
    49: 10.2,
    50: 7.7,
    51: 6.6,
    52: 5.5,
    53: 5.35,
    54: 4.044,
    55: 59.42,
    56: 39.7,
    57: 31.1,
    58: 29.6,
    59: 28.2,
    60: 31.4,
    61: 30.1,
    62: 28.8,
    63: 27.7,
    64: 23.5,
    65: 25.5,
    66: 24.5,
    67: 23.6,
    68: 22.7,
    69: 21.8,
    70: 21,
    71: 21.9,
    72: 16.2,
    73: 13.1,
    74: 11.1,
    75: 9.7,
    76: 8.5,
    77: 7.6,
    78: 6.5,
    79: 5.8,
    80: 5.02,
    81: 7.6,
    82: 6.8,
    83: 7.4,
    84: 6.8,
    85: 6,
    86: 5.3,
    87: 48.6,
    88: 38.3,
    89: 32.1,
    90: 32.1,
    91: 25.4,
    92: 24.9,
    93: 24.8,
    94: 24.5,
    95: 23.3,
    96: 23,
    97: 22.7,
    98: 20.5,
    99: 19.7,
    100: 23.8,
    101: 18.2,
    102: 17.5,
}

def compute_xdm_coeffs_erin(fn_h5, fn_atoms=None):
    with h5.File(fn_h5, 'r') as f:
        mol = IOData.from_file(f['mol'])
        moments = f['moments'][:]

    # get the (rescaled) polarizabilities
    if fn_atoms is None:
        polds = np.array([pold_map[number]*pold_unit for number in mol.numbers])
    else:
        raise NotImplementedError

    # dipole recipe
    coeffs = []
    for iatom in xrange(mol.natom):
        a = polds[iatom]
        m1 = (
            moments[  mol.natom+iatom,  mol.natom+iatom] +
            moments[2*mol.natom+iatom,2*mol.natom+iatom] +
            moments[3*mol.natom+iatom,3*mol.natom+iatom]
        )
        m2 = (
            moments[4*mol.natom+iatom,4*mol.natom+iatom] +
            moments[5*mol.natom+iatom,5*mol.natom+iatom] +
            moments[6*mol.natom+iatom,6*mol.natom+iatom] +
            moments[7*mol.natom+iatom,7*mol.natom+iatom] +
            moments[8*mol.natom+iatom,8*mol.natom+iatom]
        )
        m3 = (
            moments[ 9*mol.natom+iatom, 9*mol.natom+iatom] +
            moments[10*mol.natom+iatom,10*mol.natom+iatom] +
            moments[11*mol.natom+iatom,11*mol.natom+iatom] +
            moments[12*mol.natom+iatom,12*mol.natom+iatom] +
            moments[13*mol.natom+iatom,13*mol.natom+iatom] +
            moments[14*mol.natom+iatom,14*mol.natom+iatom] +
            moments[15*mol.natom+iatom,15*mol.natom+iatom]
        )
        denom = 2*m1*a
        c6 = m1**2*a**2/denom
        c8 = 3*(m1*m2)*a**2/denom
        c10 = (4*m1*m3 + 21.0/5.0*m2**2)*a**2/denom
        coeffs.append([c6, c8, c10])

    return np.array(coeffs)


def compute_xdm_erin(fn_h5, fn_atoms=None):
    coeffs = compute_xdm_coeffs_erin(fn_h5, fn_atoms)
    for row in coeffs:
        print '%10.5f  %10.5f  %10.5f' % tuple(row)


if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        compute_xdm_erin(sys.argv[1])
    elif len(sys.argv[1:]) == 2:
        compute_xdm_erin(sys.argv[1], sys.argv[2])
    else:
        raise NotImplementedError
