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

    # Note: these are not just the exchange holes. they also include the density itself (integrating to zero)
    exholes = []
    for iop0 in xrange(nop):
        # Construct the exchange-hole + density due with respect to iop0
        dsd = dm_full.copy()
        dsd.idot(operators[iop0])
        dsd.idot(dm_full)
        hole_dens = mol.obasis.compute_grid_density_dm(dsd, mol.grid.points)/2
        work0 = potentials[iop0]*moldens - hole_dens
        exholes.append(work0)

    moments = np.zeros((nop, nop), float)
    for iop0 in xrange(nop):
        # Integrate with other operator
        for iop1 in xrange(nop):
            moments[iop0, iop1] = mol.grid.integrate(exholes[iop0], exholes[iop1], 1.0/moldens)

    print np.diag(moments)[:4]

    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments

if __name__ == '__main__':
    main(sys.argv[1])
