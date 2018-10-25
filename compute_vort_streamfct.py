import numpy as np
import scipy.integrate as integrate  # for simpsons integration
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    timerange, kmax, nml = set_input_output_parameters(args)
    nx_half, ny_half = define_geometry(case_name, nml)
    icshift = nx_half - 1
    jcshift = ny_half -1


    # read in fields
    for it, t0 in enumerate(timerange):
        u = read_in_netcdf_fields('u', os.path.join(path_fields, str(t0) + '.nc'))
        u_roll = np.roll(np.roll(u[:, :, :], ishift, axis=0), jshift, axis=1)
        u_ = u_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        v = read_in_netcdf_fields('v', os.path.join(path_fields, str(t0) + '.nc'))
        v_roll = np.roll(np.roll(v[:, :, :], ishift, axis=0), jshift, axis=1)
        v_ = v_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        w_roll = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift, axis=1)
        w_ = w_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        # del v, w
        print 'start:', u.shape, v.shape, w.shape
        ''' FIELDS / GEOMETRY '''
        plt.figure()
        plt.subplot(131)
        plt.imshow(u_[:,:,1].T, origin='lower')
        plt.plot([icshift,icshift],[1,ny_-2],'k')
        plt.plot([1,nx_-2],[jcshift,jcshift],'k')
        plt.xlim([0,nx_])
        plt.ylim([0,ny_])
        plt.subplot(132)
        plt.imshow(v_[:,:,1].T, origin='lower')
        plt.plot([icshift, icshift], [1, ny_ - 2], 'k')
        plt.plot([1, nx_ - 2], [jcshift, jcshift], 'k')
        plt.xlim([0, nx_])
        plt.ylim([0, ny_])
        plt.subplot(133)
        plt.imshow(w_[:,:,1].T, origin='lower')
        plt.plot([icshift, icshift], [1, ny_ - 2], 'k')
        plt.plot([1, nx_ - 2], [jcshift, jcshift], 'k')
        plt.xlim([0, nx_])
        plt.ylim([0, ny_])
        plt.savefig(os.path.join(path_out, 'test_fields_t' + str(t0) + 's.png'))

        ''' VORTICITY '''
        # compute and plot vorticity in yz-cross section
        vort_yz = compute_vorticity_yz(v[ic,:,:], w[ic,:,:])
        vort_yz_ = compute_vorticity_yz(v_[icshift,:,:], w_[icshift,:,:])
        # compute and plot vorticity in xz-crosssection
        vort_xz = compute_vorticity_yz(u[:,jc,:], w[:,jc,:])
        vort_xz_ = compute_vorticity_yz(u_[:,jcshift,:], w_[:,jcshift,:])
        print('vorticity', vort_yz.shape, vort_xz.shape, nx, ny, nz)
        # >> compare the two

        fig, axes = plt.subplots(1,3, figsize=(14,4), sharey='all')
        plt.suptitle('t='+str(t0)+'s')
        ax1 = axes[0]
        a = ax1.imshow(vort_xz_.T, origin='lower')
        ax1.set_title('vort_xz')
        ax2 = axes[1]
        b = ax2.imshow(vort_yz_.T, origin='lower')
        ax2.set_title('vort_yz')
        ax3 = axes[2]
        c = ax3.imshow(vort_xz_.T-vort_yz_.T, origin='lower')
        ax3.set_title('vort_xz - vort_yz')
        plt.colorbar(a, ax=ax1)
        plt.colorbar(b, ax=ax2)
        plt.colorbar(c, ax=ax3)
        plt.savefig(os.path.join(path_out, 'vort_xz_yz_t'+str(t0)+'s.png'))
        plt.close()




        ''' STREAM FUNCTION '''
        # compute 2D streamfunction in yz-crosssection
        psi = compute_streamfunction(v[ic,:,:], w[ic,:,:], t0)

        ym = np.int(2*ny_/3)
        # compare streamfunction and streamlines
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,6))
        # plt.imshow(psi[:,:40].T, origin='lower')
        a = ax1.contourf(psi[:ym,:kmax].T)
        plt.colorbar(a, ax=ax1)
        y_arr = np.arange(0,ny)
        z_arr = np.arange(0,ny)
        ax1.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T, density=1.5, color='k')
        ax2.contourf(vort_yz[:ym,:kmax].T)
        ax2.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T)
        ax1.set_title('stream function')
        ax2.set_title('vorticity')
        plt.suptitle('t='+str(t0)+'s')
        plt.savefig(os.path.join(path_out, 'psi_t'+str(t0)+'s.png'))
        plt.close()


    return

# ----------------------------------------------------------------------
def compute_vorticity_yz(v_, w_):
    print('comp vort', v_.shape, w_.shape)
    vort_yz = np.zeros(shape=v_.shape)
    [lx, ly] = v_.shape
    for j in range(1, lx - 1):
        for k in range(1, ly-1):
            vort_yz[j, k] = (w_[j + 1, k] - w_[j - 1, k]) / (2 * dy) \
                                - (v_[j, k + 1] - v_[j, k - 1]) / (2 * dz)
    return vort_yz


def integrand(x, a, b):
    return a*x**2 + b

# def vort

def compute_streamfunction(u_, v_, t0):
    # est, errbound = integrate.quad(func, min, max)
    # funct: callable python object (function, method, class instance); also lambda-function
    # min, max: limits for integration (can us inf, -inf)
    # return:
    #   est: estimated value for integral
    #   errbound: upper bound on the error

    # testing
    # test_num_integration()

    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()

    # phi: vector potential
    # psi: stream function (solenoidal part)
    [lx, ly] = u_.shape
    print(lx, ly)

    # Integrate velocity fields to get potential and streamfunction
    # Use Simpson rule summation(function CUMSIMP)
    # for u = -dpsi/dy, v = dpsi/dx
    # psi = int(v dx) - int(u dy)
    # psi = integrate.quad()

    psi_A = np.zeros((lx, ly),dtype=np.double)
    psi_B = np.zeros((lx, ly),dtype=np.double)
    psi_simps = np.zeros((lx, ly))
    for i in range(1,lx-1):
        for j in range(1,ly-1):
            # Simpson Rule
            psi_A[i,j] = 2*dx/6 * ( v_[i+1,j] + 4*v_[i,j] + v_[i-1,j] )
            psi_B[i,j] = 2*dy/6 * ( u_[i,j+1] + 4*u_[i,j] + u_[i,j-1] )
    psi_simps = psi_A - psi_B

    km = 30
    ym = np.int(ny_/2)
    fig, axes = plt.subplots(1,2, figsize=(10,4))
    ax1 = plt.subplot(121)
    ax1.set_title('u')
    plt.contourf(u_[:ym,:km].T, alpha=0.7)
    plt.colorbar()
    plt.contour(u_[:ym,:km].T, levels=[0], colors='k', linewidth=5)
    CS = ax1.contour(-psi_B[:ym,:km].T, linewidth=3)
    ax1.clabel(CS)
    ax2 = plt.subplot(122)
    ax2.set_title('psi_B = int dy u')
    print('maxs', np.amax(u_[:ym,:km]), np.amax(psi_B[:ym,:km]))
    ax2.contour(-psi_B[:ym,:km].T)
    ax2.contour(u_[:ym,:km].T, linestyles='--')
    fig.savefig(os.path.join(path_out, 'test_u_t'+str(t0)+'s.png'))
    plt.close()


    return psi_simps


def test_num_integration():
    x = np.arange(1, 10)
    a = 2
    f = a * x
    I = integrate.simps(f, x)
    print x
    print f
    print I
    x = np.arange(0, 10)
    f = 2 * np.ones(10)
    I = integrate.simps(f, x)
    print f
    print I
    x = np.arange(0, 10)
    y = np.arange(0, 5)
    X, Y = np.meshgrid(x, y, indexing='ij')  # if indexing not defined, output will be 'xy'-indexed
    # 'ij'-indexing: X.shape = Y.shape = [len(x), len(y)]
    # 'xy'-indexing: X.shape = Y.shape = [len(y), len(x)]
    # X[:,j] = x, Y[i,:] = y
    a = 1
    b = 2
    f = a * X + b * Y
    print X
    print Y
    print f

    i = 0
    j = 0
    I = integrate.simps(f[:, j], X[:, j])
    print I
    I = integrate.simps(f[:, j], x, axis=0)
    print I
    I = integrate.simps(f, x, axis=0)
    print I
    I2 = integrate.simps(f, y, axis=1)
    print I2


    a = 2
    b = 1
    # integrand(x, 1, 0)
    psi = integrate.quad(integrand, 0, 1, args=(a,b))
    print('psi test', psi)


    return

# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('--- set input parameters ---')
    global path_in, path_out, path_stats, path_fields
    if args.path:
        path = args.path
    else:
        # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
        #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th10K/'
    path_in = os.path.join(path, 'fields_cp_rim')
    path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'figs_vorticity')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx * dy * dz

    global nt
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)

    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 50

    return timerange, kmax, nml



def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global nx_, ny_
    global ic, jc, shift, ishift, jshift
    if case_name == 'ColdPoolDry_triple_3D':
        flag = 'triple'
        # d = np.int(np.round(ny / 2))
        d = np.int(np.round((ny + gw) / 2))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        ic = np.int(np.round(a / 2))
        jc = np.int(np.round(d / 2))
        shift = 60
        nx_half = irstar + shift
        ny_half = irstar + shift
        ishift = np.max(nx_half - ic, 0)
        jshift = np.max(ny_half - jc, 0)
        nx_ = 2 * nx_half
        ny_ = 2 * ny_half
    elif case_name == 'ColdPoolDry_double_3D':
        flag = 'double'
        rstar = 5000.0
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(nx / 3)
        jc1 = np.int(ny / 2)
        # ic2 = ic1 + isep
        # jc2 = jc1 + jsep
        ic = ic1
        jc = jc1
        shift = 40
        nx_half = irstar + shift
        ny_half = irstar + shift
        ishift = np.max(nx_half - ic, 0)
        jshift = np.max(ny_half - jc, 0)
        nx_ = 2 * nx_half
        ny_ = 2 * ny_half

    print('rstar: ' + str(rstar), irstar)
    print('ic,jc,id,jd', ic, jc, nx_half, ny_half)
    print('nx_,ny_', nx_, ny_)
    print('shift, ishift, jshift', shift, ishift, jshift)
    return nx_half, ny_half


# ----------------------------------------------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()