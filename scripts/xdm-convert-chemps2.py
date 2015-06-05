#!/usr/bin/env python

import sys
import h5py as h5, numpy as np
from horton import *


def load_molden_irreps(fn_molden):
    irreps = []
    with open(fn_molden) as f:
        for line in f:
            if line.startswith(' Sym='):
                irreps.append(line[5:].strip())
            if line == ' Spin= Beta\n':
                del irreps[-1]
                break
    return np.array(irreps)


def convert_chemps2(fn_2dm, fn_hf, fn_out):
    # Load the 2RDM
    with h5.File(fn_2dm, 'r') as f:
        energy = f['energy'][()]
        twodm_mo = f['twodm'][:]/2
    dt = np.einsum('abab', twodm_mo)
    #print 'dt', dt
    nel = 0.5 * ( 1.0 + np.sqrt(1.0 + 8*dt) )
    #print 'nel', nel
    onedm_mo = np.einsum('acbc->ab', twodm_mo)*2/(nel-1)
    #print np.trace(onedm_mo)

    # Load the HF orbitals to transform the 2RDM back to the orbital basis
    mol = IOData.from_file(fn_hf)
    mol.grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'ultrafine', mode='keep', random_rotate=False)

    # Load the irreps from the molden file
    irreps_molden = load_molden_irreps(fn_hf)
    print irreps_molden
    print len(irreps_molden)
    irreps_order_chemps2 = ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']

    # Construct the Chemps2 order of the orbitals:
    order = np.zeros(len(irreps_molden), int)
    counter = 0
    for irrep in irreps_order_chemps2:
        for i in xrange(len(irreps_molden)):
            if irreps_molden[i] == irrep:
                order[counter] = i
                counter += 1
    print order
    print irreps_molden[order]

    # Reverse transformation: from mo to ao
    exp = mol.exp_alpha
    print exp.coeffs.shape
    tf = exp.coeffs[:,order].T

    twodm_ao = np.tensordot(twodm_mo, tf, axes=([0],[0]))
    twodm_ao = np.tensordot(twodm_ao, tf, axes=([0],[0]))
    twodm_ao = np.tensordot(twodm_ao, tf, axes=([0],[0]))
    twodm_ao = np.tensordot(twodm_ao, tf, axes=([0],[0]))

    onedm_ao = np.tensordot(onedm_mo, tf, axes=([0],[0]))
    onedm_ao = np.tensordot(onedm_ao, tf, axes=([0],[0]))

    # conversion and final check
    lf = DenseLinalgFactory(mol.obasis.nbasis)
    olp = mol.obasis.compute_overlap(mol.lf)
    onedm_ao_horton = lf.create_two_index()
    onedm_ao_horton._array[:] = onedm_ao
    print olp.contract_two('ab,ab', onedm_ao_horton)

    twodm_ao_horton = lf.create_four_index()
    twodm_ao_horton._array[:] = twodm_ao

    tmp0 = twodm_ao_horton.contract_two_to_two('abcd,bd->ac', olp)
    print olp.contract_two('ab,ab', tmp0)

    # Write all to a HDF file
    with h5.File(fn_out, 'w') as f:
        f.create_group('mol')
        mol.to_file(f['mol'])
        del f['mol']['exp_alpha']
        f['mol'].create_group('dm_full')
        dump_h5(f['mol/dm_full'], onedm_ao_horton)
        f['mol'].create_group('twodm_full')
        dump_h5(f['mol/twodm_full'], twodm_ao_horton)

if __name__ ==  '__main__':
    convert_chemps2(sys.argv[1], sys.argv[2], sys.argv[3])
