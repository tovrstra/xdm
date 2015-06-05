#!/usr/bin/env python

import sys
import h5py as h5
from horton import *

if __name__ == '__main__':
    mol = IOData.from_file(sys.argv[1])
    mol.grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'ultrafine', mode='keep', random_rotate=False)
    with h5.File(sys.argv[2], 'w') as f:
        f.create_group('mol')
        mol.to_file(f['mol'])
