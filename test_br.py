#!/usr/bin/env python

from cext import solve_brhole_isotropic
import numpy as np, h5py as h5, matplotlib.pyplot as pt
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
        A = brhole[:,0]
        a = brhole[:,1]
        b = brhole[:,2]

        begin = 0
        for iatom in xrange(mol.natom):
            p = mol.grid.points - mol.coordinates[iatom]
            r = (p[:,0]*p[:,0] + p[:,1]*p[:,1] + p[:,2]*p[:,2])**0.5
            rb = r - b
            #rb[b==0] = 0.0
            #rb[:] = 0.0
            w = potentials['mom_00_atom_%03i' % iatom]

            atgrid = mol.grid.subgrids[iatom]
            end = begin + atgrid.size

            for imom in 3, 2, 1:
                xd = r**imom - rb**imom
                #xd[b==0] = 0.0
                moments[iatom, imom-1] = mol.grid.integrate(w, xd, xd, mgga[:,0])

                pt.clf()
                int_sav = atgrid.get_spherical_average(w[begin:end], xd[begin:end], xd[begin:end], mgga[begin:end,0])
                pt.plot(atgrid.rgrid.radii, np.cumsum(atgrid.rgrid.weights*int_sav), label='int')
                #pt.legend(loc=0)
                pt.xlim(0, 10)
                pt.savefig('integrand_%03i_%1i.png' % (iatom, imom))

            r_sav = atgrid.get_spherical_average(r[begin:end])
            b_sav = atgrid.get_spherical_average(b[begin:end])
            rb_sav = atgrid.get_spherical_average(rb[begin:end])
            a_sav = atgrid.get_spherical_average(a[begin:end])
            A_sav = atgrid.get_spherical_average(A[begin:end])
            rho_sav = atgrid.get_spherical_average(mgga[begin:end,0])
            rgrad_sav = atgrid.get_spherical_average((mgga[begin:end,1]**2+mgga[begin:end,2]**2+mgga[begin:end,3]**2)/mgga[begin:end,0])
            lapl_sav = atgrid.get_spherical_average(mgga[begin:end,4])
            tau_sav = atgrid.get_spherical_average(mgga[begin:end,5])

            pt.clf()
            pt.plot(atgrid.rgrid.radii, r_sav, label='r')
            pt.plot(atgrid.rgrid.radii, b_sav, label='b')
            pt.plot(atgrid.rgrid.radii, rb_sav, label='rb')
            pt.legend(loc=0)
            pt.xlim(0,10)
            pt.ylim(-1,10)
            pt.savefig('radii_%03i.png' % iatom)

            pt.clf()
            pt.semilogy(atgrid.rgrid.radii, rho_sav, 'b-',label='rho')
            pt.semilogy(atgrid.rgrid.radii, rgrad_sav, 'g-',label='rgrad')
            pt.semilogy(atgrid.rgrid.radii, lapl_sav, 'r-', label='lapl')
            pt.semilogy(atgrid.rgrid.radii, -lapl_sav, 'r--', label='lapl')
            pt.semilogy(atgrid.rgrid.radii, tau_sav, 'k-', label='tau')
            pt.legend(loc=0)
            pt.xlim(0,10)
            pt.ylim(1e-12,1e3)
            pt.savefig('rho_%03i.png' % iatom)

            pt.clf()
            pt.semilogy(atgrid.rgrid.radii, A_sav, label='A')
            pt.semilogy(atgrid.rgrid.radii, a_sav, label='a')
            pt.axhline(1/np.pi, color='b', ls='--')
            pt.axhline(2.0, color='g', ls='--')
            pt.legend(loc=0)
            pt.xlim(0,10)
            pt.savefig('slater_%03i.png' % iatom)

            begin = end


    if restricted:
        moments *= 2

    print moments[0]
    return

    print moments
    with h5.File(fn_work) as f:
        if 'moments' in f:
            del f['moments']
        f['moments'] = moments


if __name__ == '__main__':
    main(sys.argv[1])
