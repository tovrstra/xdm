#!/usr/bin/env python

from horton import *
import numpy as np, h5py as h5
import sys

log.set_level(log.silent)

def main(fn_work):
    with h5.File(fn_work, 'r') as f:
        mol = IOData.from_file(f['mol'])
        mol.lf = DenseLinalgFactory(mol.obasis.nbasis)
        operators = load_h5(f['aim']['operators'])
        potentials = load_h5(f['aim']['potentials'])

    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, mol.grid.points)

    nop = len(operators)
    natom = mol.natom
    npure = nop/natom

    operators = [value for key, value in sorted(operators.iteritems())]
    potentials = [value for key, value in sorted(potentials.iteritems())]

    moments = np.zeros((natom, 3), float)
    for iatom in xrange(natom):
        bary_centers = []
        pots = []
        for ipure in xrange(1, 4):
            iop = iatom + ipure*natom
            # Expectation value of the X-hole barycenter
            dsd = dm_full.copy()
            dsd.idot(operators[iop])
            dsd.idot(dm_full)
            bary_center = mol.obasis.compute_grid_density_dm(dsd, mol.grid.points)/2/moldens
            bary_centers.append(bary_center)
            pots.append(potentials[iop])


        cmag = 0.0
        rmag = 0.0
        for icart in 0, 1, 2:
            cmag += bary_centers[icart]**2
            rmag += pots[icart]**2
        cmag **= 0.5
        rmag **= 0.5

        for imom in 1, 2, 3:
            xd = rmag**imom - cmag**imom
            moments[iatom, imom-1] = mol.grid.integrate(xd, xd, moldens)

    print moments[0]

    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments

if __name__ == '__main__':
    main(sys.argv[1])
