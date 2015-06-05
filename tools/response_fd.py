#!/usr/bin/env python

from horton import *
import numpy as np
import sys


def main(fn_mol):
    mol = IOData.from_file(fn_mol)
    mol.obasis = get_gobasis(mol.coordinates, mol.numbers, 'aug-cc-pvtz')
    mol.lf = DenseLinalgFactory(mol.obasis.nbasis)

    # Compute Gaussian integrals
    mol.olp = mol.obasis.compute_overlap(mol.lf)
    mol.kin = mol.obasis.compute_kinetic(mol.lf)
    mol.na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    mol.er = mol.obasis.compute_electron_repulsion(mol.lf)

    # Setup integration grids with default settings
    mol.grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='keep')

    # Ground state calculation
    dm_alpha = do_scf(mol)

    # MBIS partitioning
    moldens = 2*mol.obasis.compute_grid_density_dm(dm_alpha, mol.grid.points)
    wpart = MBISWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, mol.grid, moldens, lmax=1)
    wpart.do_charges()
    print wpart['charges']

    # Construct the distributed multipole operators
    operators, potentials = get_dmo(mol, wpart)
    nop = len(operators)

    # Do finite differences to construct the response matrix
    response_matrix = np.zeros((nop, nop), float)
    eps = 1e-3
    for iop0 in xrange(nop):
        # perturbation in one direction
        pert = operators[iop0].copy()
        pert.iscale(eps)
        dm_alpha_plus = dm_alpha.copy()
        do_scf(mol, dm_alpha_plus, pert)
        # perturbation in other direction
        pert.iscale(-1)
        dm_alpha_min = dm_alpha.copy()
        do_scf(mol, dm_alpha_min, pert)
        # the response density matrix
        dm_resp = dm_alpha_plus.copy()
        dm_resp.iadd(dm_alpha_min, -1)
        dm_resp.iscale(1.0/eps)
        # compute matrix elements
        for iop1 in xrange(nop):
            response_matrix[iop0,iop1] = dm_resp.contract_two('ab,ab', operators[iop1])
    # symmetrize
    response_matrix = (response_matrix + response_matrix.T)/2

    # write stuff out to a HDF5 file. Keep as much as possible
    results = {}
    npure = get_npure_cumul(wpart.lmax)
    for key, value in mol.__dict__.iteritems():
        if key == 'lf':
            continue
        if key.startswith('_'):
            key = key[1:]
        results[key] = value
    results['operators'] = {}
    results['potentials'] = {}
    for ipure in xrange(npure):
        for iatom in xrange(mol.natom):
            results['operators']['mom_%02i_atom_%03i' % (ipure, iatom)] = operators.pop(0)
            results['potentials']['mom_%02i_atom_%03i' % (ipure, iatom)] = potentials.pop(0)
    dump_h5('foo.h5', results)

def get_dmo(mol, wpart):
    wpart.do_partitioning()
    npure = get_npure_cumul(wpart.lmax)
    for index in xrange(wpart.natom):
        if log.do_medium:
            log('Computing overlap matrices for atom %i.' % index)

        # Prepare solid harmonics on grids.
        grid = wpart.get_grid(index)
        if wpart.lmax > 0:
            work = np.zeros((grid.size, npure-1), float)
            work[:,0] = grid.points[:,2] - wpart.coordinates[index,2]
            work[:,1] = grid.points[:,0] - wpart.coordinates[index,0]
            work[:,2] = grid.points[:,1] - wpart.coordinates[index,1]
            if wpart.lmax > 1:
                fill_pure_polynomials(work, wpart.lmax)
        at_weights = wpart.cache.load('at_weights', index)

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

        wpart.cache.dump(('overlap_operators', index), overlap_operators)
        wpart.cache.dump(('potential_functions', index), potential_functions)

    # Correct the s-type overlap operators such that the sum is exactly
    # equal to the total overlap.
    error_overlap = mol.lf.create_two_index()
    for index in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', index)['olp_00000']
        error_overlap.iadd(atom_overlap)
    error_overlap.iadd(mol.obasis.compute_overlap(mol.lf), -1)
    error_overlap.iscale(1.0/wpart.natom)
    for index in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', index)['olp_00000']
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


def do_scf(mol, dm_alpha0=None, pert=None):
    # Construction of Hamiltonian
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    libxc_term = RLibXCHybridGGA('xc_b3lyp')
    terms = [
        RTwoIndexTerm(mol.kin, 'kin'),
        RDirectTerm(mol.er, 'hartree'),
        RGridGroup(mol.obasis, mol.grid, [
            libxc_term,
        ]),
        RExchangeTerm(mol.er, 'x_hf', libxc_term.get_exx_fraction()),
        RTwoIndexTerm(mol.na, 'ne'),
    ]
    if pert is not None:
        terms.append(RTwoIndexTerm(pert, 'pert'))
    ham = REffHam(terms, external)

    # Decide how to occupy the orbitals
    occ_model = AufbauOccModel(mol.numbers.sum()/2)

    # Construct the initial density matrix
    if dm_alpha0 is not None:
        dm_alpha = dm_alpha0
    else:
        exp_alpha0 = mol.lf.create_expansion()
        guess_core_hamiltonian(mol.olp, mol.kin, mol.na, exp_alpha0)
        occ_model.assign(exp_alpha0)
        dm_alpha = exp_alpha0.to_dm()

    # Optimal damping SCF cycles
    scf_solver = ODASCFSolver(1e-6)
    scf_solver(ham, mol.lf, mol.olp, occ_model, dm_alpha)

    # Return the optimized density matrix
    return dm_alpha


if __name__ == '__main__':
    main(sys.argv[1])
