#!/usr/bin/env python

from horton import *
import numpy as np, h5py as h5
import sys


def get_volume(spline, pro=False):
    rgrid = RadialGrid(spline.rtransform)
    rho = spline.y
    if not pro:
        rho /= np.sqrt(4*np.pi)
    r = rgrid.radii
    return rgrid.integrate(rho, r**3)/rgrid.integrate(rho)


def main(fn_work, scheme, fn_atoms=None):
    with h5.File(fn_work, 'r') as f:
        mol = IOData.from_file(f['mol'])
    mol.lf = DenseLinalgFactory(mol.obasis.nbasis)

    moldens = mol.obasis.compute_grid_density_dm(mol.get_dm_full(), mol.grid.points)

    padb = ProAtomDB.from_file(fn_atoms)
    if scheme == 'mbis':
        wpart = MBISWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, mol.grid, moldens, lmax=3)
    elif scheme == 'hi':
        wpart = HirshfeldIWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, mol.grid, moldens, padb, lmax=3)
    elif scheme == 'h':
        wpart = HirshfeldWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, mol.grid, moldens, padb, lmax=3)
    else:
        raise Foo

    wpart.do_charges()
    wpart.do_moments()

    # compute the ratio of the expectation value of r^3
    volume_ratios = []
    for iatom in xrange(mol.natom):
        spline = padb.get_spline(mol.numbers[iatom])
        pro_volume = get_volume(spline, pro=True)
        assert abs(wpart['radial_moments'][iatom, 0] - wpart['populations'][iatom]) < 1e-10
        aim_volume = wpart['radial_moments'][iatom, 3]/wpart['radial_moments'][iatom, 0]
        volume_ratios.append(aim_volume/pro_volume)
    volume_ratios = np.array(volume_ratios)

    # Construct the distributed multipole operators
    operators, potentials = get_dmo(mol, wpart)
    nop = len(operators)

    # write stuff out to a HDF5 file. Keep as much as possible
    results = {}
    npure = get_npure_cumul(wpart.lmax)
    results['volume_ratios'] = volume_ratios
    results['charges'] = wpart['charges']
    results['radial_moments'] = wpart['radial_moments']
    results['operators'] = {}
    results['potentials'] = {}
    for ipure in xrange(npure):
        for iatom in xrange(mol.natom):
            results['operators']['mom_%02i_atom_%03i' % (ipure, iatom)] = operators.pop(0)
            results['potentials']['mom_%02i_atom_%03i' % (ipure, iatom)] = potentials.pop(0)

    with h5.File(fn_work) as f:
        if 'aim' in f:
            del f['aim']
        f.create_group('aim')
        dump_h5(f['aim'], results)


def get_dmo(mol, wpart):
    wpart.do_partitioning()
    npure = get_npure_cumul(wpart.lmax)

    # Build the AIM weight function for the total grid
    promol = mol.grid.zeros()
    all_at_weights = []
    for iatom in xrange(wpart.natom):
        proatom = mol.grid.zeros()
        wpart.eval_proatom(iatom, proatom, mol.grid)
        all_at_weights.append(proatom)
        promol += proatom
    for iatom in xrange(wpart.natom):
        all_at_weights[iatom] /= promol
    del promol

    # Build operators
    for iatom in xrange(wpart.natom):
        if log.do_medium:
            log('Computing overlap matrices for atom %i.' % iatom)

        # Prepare solid harmonics on grids.
        grid = mol.grid
        if wpart.lmax > 0:
            work = np.zeros((grid.size, npure-1), float)
            work[:,0] = grid.points[:,2] - wpart.coordinates[iatom,2]
            work[:,1] = grid.points[:,0] - wpart.coordinates[iatom,0]
            work[:,2] = grid.points[:,1] - wpart.coordinates[iatom,1]
            if wpart.lmax > 1:
                fill_pure_polynomials(work, wpart.lmax)

        at_weights = all_at_weights[iatom]

        # Convert the weight functions to AIM overlap operators.
        counter = 0
        overlap_operators = {}
        potential_functions = {}
        for l in xrange(wpart.lmax+1):
            for m in xrange(-l, l+1):
                op = mol.lf.create_two_index()
                if counter > 0:
                    tmp = at_weights*work[:,counter-1]
                else:
                    tmp = at_weights
                mol.obasis.compute_grid_density_fock(grid.points, grid.weights, tmp, op)
                overlap_operators['olp_%05i' % counter] = op
                potential_functions['pot_%05i' % counter] = tmp
                counter += 1

        wpart.cache.dump(('overlap_operators', iatom), overlap_operators)
        wpart.cache.dump(('potential_functions', iatom), potential_functions)

    # Correct the s-type overlap operators such that the sum is exactly
    # equal to the total overlap.
    error_overlap = mol.lf.create_two_index()
    for iatom in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', iatom)['olp_00000']
        error_overlap.iadd(atom_overlap)
    error_overlap.iadd(mol.obasis.compute_overlap(mol.lf), -1)
    error_overlap.iscale(1.0/wpart.natom)
    for iatom in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', iatom)['olp_00000']
        atom_overlap.iadd(error_overlap, -1)

    # Construct the list of operators with a logical ordering
    #   * outer loop: s, pz, px, py, ...
    #   * inner loop: atoms
    operators = []
    potentials = []
    for ipure in xrange(npure):
        for iatom in xrange(wpart.natom):
            operators.append(
                wpart.cache.load('overlap_operators', iatom)['olp_%05i' % ipure]
            )
            potentials.append(
                wpart.cache.load('potential_functions', iatom)['pot_%05i' % ipure],
            )
    return operators, potentials


if __name__ == '__main__':
    fn_atoms = None
    if len(sys.argv) == 4:
        fn_atoms = sys.argv[3]
    main(sys.argv[1], sys.argv[2], fn_atoms)
