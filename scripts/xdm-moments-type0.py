#!/usr/bin/env python

from horton import *
import numpy as np, h5py as h5
import sys

log.set_level(log.silent)

def main(fn_work):
    with h5.File(fn_work, 'r') as f:
        mol = IOData.from_file(f['mol'])
        mol.lf = DenseLinalgFactory(mol.obasis.nbasis)
        potentials = load_h5(f['aim']['potentials'])

    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, mol.grid.points)

    field_x = mol.lf.create_two_index()
    field_y = mol.lf.create_two_index()
    field_z = mol.lf.create_two_index()
    mol.obasis.compute_grid_density_fock(mol.grid.points, mol.grid.weights, mol.grid.points[:,0], field_x)
    mol.obasis.compute_grid_density_fock(mol.grid.points, mol.grid.weights, mol.grid.points[:,1], field_y)
    mol.obasis.compute_grid_density_fock(mol.grid.points, mol.grid.weights, mol.grid.points[:,2], field_z)

    # Compute x-hole dipole components
    xds = np.zeros(mol.grid.points.shape)
    for icart, op in enumerate([field_x, field_y, field_z]):
        # Expectation value of the X-hole barycenter
        dsd = dm_full.copy()
        dsd.idot(op)
        dsd.idot(dm_full)
        bary_center = mol.obasis.compute_grid_density_dm(dsd, mol.grid.points)/2/moldens
        xds[:,icart] = mol.grid.points[:,icart] - bary_center

    # Compute the magnitude of the x-hole dipole
    xdmag = np.sqrt(xds[:,0]*xds[:,0] + xds[:,1]*xds[:,1] + xds[:,2]*xds[:,2])

    moments = np.zeros((mol.natom, 3), float)
    for iatom in xrange(mol.natom):
        r = mol.grid.points - mol.coordinates[iatom]
        rmag = np.sqrt(r[:,0]*r[:,0] + r[:,1]*r[:,1] + r[:,2]*r[:,2])
        w = potentials['mom_00_atom_%03i' % iatom]
        for imom in 1, 2, 3:
            tmp = rmag**imom - (rmag - xdmag)**imom
            moments[iatom, imom-1] = mol.grid.integrate(w, tmp, tmp, moldens)

    print moments[0]

    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments

if __name__ == '__main__':
    main(sys.argv[1])
