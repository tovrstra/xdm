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

    onedm_full = mol.dm_full
    twodm_full = mol.twodm_full
    moldens = mol.obasis.compute_grid_density_dm(onedm_full, mol.grid.points)

    nop = len(operators)
    operators = [value for key, value in sorted(operators.iteritems())]
    potentials = [value for key, value in sorted(potentials.iteritems())]

    # expecation values of all the operators
    expvals = [onedm_full.contract_two('ab,ab', operator) for operator in operators]
    print expvals[:4]

    moments = np.zeros((nop, nop), float)
    for iop0 in xrange(nop):
        work0 = potentials[iop0]*moldens
        tmp0 = twodm_full.contract_two_to_two('abcd,bd->ac', operators[iop0])

        # Integrate with other operator
        for iop1 in xrange(nop):
            moments[iop0, iop1] = sum([
                2*tmp0.contract_two('ab,ab', operators[iop1]),
                +mol.grid.integrate(work0, potentials[iop1]),
                -expvals[iop0]*expvals[iop1],
            ])

    print np.diag(moments)[:4]

    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments

if __name__ == '__main__':
    main(sys.argv[1])
