#!/usr/bin/env python

from cext import solve_brhole_isotropic
import numpy as np, h5py as h5
import sys

from horton import *

log.set_level(log.silent)


def main(fn_work):
    with h5.File(fn_work, 'r') as f:
        mol = IOData.from_file(f['mol'])
        mol.lf = DenseLinalgFactory(mol.obasis.nbasis)
        operators = load_h5(f['aim']['operators'])
        potentials = load_h5(f['aim']['potentials'])

    if hasattr(mol, 'exp_beta'):
        exps = [mol.exp_alpha, mol.exp_beta]
        restricted = False
    else:
        exps = [mol.exp_alpha]
        restricted = True

    moments = np.zeros((mol.natom, 3))
    for exp in exps:
        if sum(exp.occupations) == 0.0:
            print 'Skipping empty spin channel'
            continue
        dm = exp.to_dm()
        mgga = mol.obasis.compute_grid_mgga_dm(dm, mol.grid.points)

        brhole = np.array([solve_brhole_isotropic(row, threshold=1e-6) for row in mgga])
        #brhole = solve_brhole_isotropic_ar ray(mgga)
        b = brhole[:,2]

        for iatom in xrange(mol.natom):
            p = mol.grid.points - mol.coordinates[iatom]
            r = (p[:,0]*p[:,0] + p[:,1]*p[:,1] + p[:,2]*p[:,2])**0.5
            rb = r - b
            #rb[b==0] = 0.0
            #rb[:] = 0.0
            w = potentials['mom_00_atom_%03i' % iatom]
            #print mgga[::1000,0]*2
            for imom in 1, 2, 3:
                xd = r**imom - rb**imom
                #xd[b==0.0] = 0.0
                moments[iatom, imom-1] += mol.grid.integrate(w, xd, xd, mgga[:,0])

    if restricted:
        moments *= 2

    print moments[0]
    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments


if __name__ == '__main__':
    main(sys.argv[1])
