#!/usr/bin/env python

import sys, argparse
from horton import *
import h5py as h5, numpy as np

polqs = {
    1: 15.0,
    2: 2.44,
    10: 6.42,
    18: 48.2,
    36: 78.8,
}

polos = {
    1: 131.25,
    2: 10.0,
    10: 30.4,
    18: 450.0,
    36: 681.0,
}


class XDMMolecule(object):
    def __init__(self, fn_h5):
        with h5.File(fn_h5, 'r') as f:
            mol = IOData.from_file(f['mol'])
            self.natom = mol.natom
            self.numbers = mol.numbers
            if f['moments'].shape == (self.natom, 3):
                self.moments = f['moments'][:]
            else:
                moments = f['moments'][:]
                self.moments = np.zeros((self.natom, 3), float)
                for iatom in xrange(self.natom):
                    m1 = (
                        moments[  self.natom+iatom,  self.natom+iatom] +
                        moments[2*self.natom+iatom,2*self.natom+iatom] +
                        moments[3*self.natom+iatom,3*self.natom+iatom]
                    )
                    m2 = (
                        moments[4*self.natom+iatom,4*self.natom+iatom] +
                        moments[5*self.natom+iatom,5*self.natom+iatom] +
                        moments[6*self.natom+iatom,6*self.natom+iatom] +
                        moments[7*self.natom+iatom,7*self.natom+iatom] +
                        moments[8*self.natom+iatom,8*self.natom+iatom]
                    )
                    m3 = (
                        moments[ 9*self.natom+iatom, 9*self.natom+iatom] +
                        moments[10*self.natom+iatom,10*self.natom+iatom] +
                        moments[11*self.natom+iatom,11*self.natom+iatom] +
                        moments[12*self.natom+iatom,12*self.natom+iatom] +
                        moments[13*self.natom+iatom,13*self.natom+iatom] +
                        moments[14*self.natom+iatom,14*self.natom+iatom] +
                        moments[15*self.natom+iatom,15*self.natom+iatom]
                    )
                    self.moments[iatom] = m1, m2, m3

            # get the (rescaled) polarizabilities
            self.pols = []
            for iatom in xrange(mol.natom):
                number = mol.numbers[iatom]
                self.pols.append([
                    periodic[number].pold*f['aim/volume_ratios'][iatom],
                    polqs.get(number, 0.0),
                    polos.get(number, 0.0)
                ])
            self.pols = np.array(self.pols)


def compute_xdm_coeffs_erin(ma, mb, pa, pb, dqo=False):
    if dqo:
        fac23 = 2.0/3.0
        fac25 = 2.0/5.0
        fac27 = 2.0/7.0
        denom_dd    = fac23*ma[0]/pa[0] + fac23*mb[0]/pb[0]
        denom_qd_ab = fac25*ma[1]/pa[1] + fac23*mb[0]/pb[0]
        denom_qd_ba = fac23*ma[0]/pa[0] + fac25*mb[1]/pb[1]
        denom_qq    = fac25*ma[1]/pa[1] + fac25*mb[1]/pb[1]
        denom_od_ab = fac27*ma[2]/pa[2] + fac23*mb[0]/pb[0]
        denom_od_ba = fac23*ma[0]/pa[0] + fac27*mb[2]/pb[2]
        c6 = fac23*ma[0]*mb[0]/denom_dd
        c8 = ma[1]*mb[0]/denom_qd_ab + mb[1]*ma[0]/denom_qd_ba
        c10 = 4.0/3.0*(ma[2]*mb[0]/denom_od_ab + mb[2]*ma[0]/denom_od_ba) + 2.8*ma[1]*mb[1]/denom_qq
        print c6, c8, c10
    else:
        factor = pa[0]*pb[0]/(ma[0]*pb[0]+mb[0]*pa[0])
        c6 = factor*ma[0]*mb[0]
        c8 = factor*1.5*(ma[1]*mb[0] + mb[1]*ma[0])
        c10 = factor*(2*(ma[2]*mb[0] + mb[2]*ma[0]) + 4.2*ma[1]*mb[1])
    return c6, c8, c10


def compute_xdm_erin(fn1_h5, fn2_h5, fn_out, dqo=False):
    mol1 = XDMMolecule(fn1_h5)
    mol2 = XDMMolecule(fn2_h5)
    coeffs = np.zeros((mol1.natom, mol2.natom, 3), float)
    for iatom1 in xrange(mol1.natom):
        for iatom2 in xrange(mol2.natom):
            coeffs[iatom1,iatom2] = \
                compute_xdm_coeffs_erin(mol1.moments[iatom1], mol2.moments[iatom2], mol1.pols[iatom1], mol2.pols[iatom2], dqo)
    with h5.File(fn_out, 'w') as f:
        f['c6_coeffs'] = coeffs[:,:,0]
        f['c8_coeffs'] = coeffs[:,:,1]
        f['c10_coeffs'] = coeffs[:,:,2]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='xdm_erin.py',
        description='Compute C6, C8 and C10 coeffs using Erin\'s XDM')

    parser.add_argument('input1',
        help='The first HDF5 file with moments')
    parser.add_argument('input2',
        help='The second HDF5 file with moments')
    parser.add_argument('output',
        help='The output HDF5 file')
    parser.add_argument('--dqo', default=False, action='store_true',
        help='Use higher multipole polarizabilities.')

    args = parser.parse_args()
    compute_xdm_erin(args.input1, args.input2, args.output, args.dqo)
