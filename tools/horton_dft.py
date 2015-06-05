#!/usr/bin/env python

from horton import *
import numpy as np, h5py as h5
import sys


def main(fn_mol, fn_out):
    mol = IOData.from_file(fn_mol)
    setup_molecule(mol)
    dm_full = do_scf(mol)
    dm_full.iscale(2)

    with h5.File(fn_out, 'w') as f:
        f.create_group('mol')
        mol.to_file(f['mol'])
        f['mol'].create_group('dm_full')
        dump_h5(f['mol/dm_full'], dm_full)



def setup_molecule(mol):
    mol.obasis = get_gobasis(mol.coordinates, mol.numbers, 'aug-cc-pvtz')
    mol.lf = DenseLinalgFactory(mol.obasis.nbasis)

    # Compute Gaussian integrals
    mol.olp = mol.obasis.compute_overlap(mol.lf)
    mol.kin = mol.obasis.compute_kinetic(mol.lf)
    mol.na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    mol.er = mol.obasis.compute_electron_repulsion(mol.lf)

    # Setup integration grids with default settings
    mol.grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'ultrafine', mode='keep', random_rotate=False)


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
    main(sys.argv[1], sys.argv[2])
