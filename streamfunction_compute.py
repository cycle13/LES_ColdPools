import numpy as np
import scipy.integrate as integrate  # for simpsons integration
import scipy.integrate as integrate  # for simpsons integration
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys
import time

'''
COMPUTE STREAMFUNCTION

'''


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    timerange, kmax, nml = set_input_output_parameters(args)

    ''' TEST FUNCTIONS '''
    # streamfunction_circular_flow()

    ''' STREAM FUNCTION '''
    # for it, t0 in enumerate(timerange):
    #     print('--- time: t=' + str(t0) + 's ---')
    #
    #     # compute 2D streamfunction in yz-crosssection
    #     psi_1 = compute_streamfunction_simps(v[ic,:,:], w[ic,:,:], t0)
    #     psi_2 = compute_streamfunction_simps2(v[ic,:,:], w[ic,:,:], t0)
    #     psi_3 = compute_streamfunction_simps3(v[ic,:,:], w[ic,:,:], t0)
    #     psi_4 = compute_streamfunction_simps4(v[ic,:,:], w[ic,:,:], rho0, t0)
    #     psi_5 = compute_streamfunction_simps5(v[ic,:,:], w[ic,:,:], rho0, 0, 0)
    #
    #     # compare two functions
    #     ymax = 120
    #     y = np.arange(0,ymax)
    #     z = np.arange(kmax)
    #     Y, Z = np.meshgrid(y, z,indexing='ij')
    #     print 'psi', psi_1.shape, psi_2.shape, y.shape, z.shape, Y.shape, Z.shape
    #     fig, axes = plt.subplots(2,3, figsize=(16,8))
    #     ax1 = axes[0,0]
    #     ax1.set_title('psi_1')
    #     max = np.maximum(np.amax(psi_1[:ymax, :kmax]),-np.amin(psi_1[:ymax, :kmax]))
    #     a = ax1.contourf(Y,Z,psi_1[:ymax, :kmax], levels=np.linspace(-max,max,1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax2 = axes[0,1]
    #     ax2.set_title('psi_2')
    #     max = np.maximum(np.amax(psi_2[:ymax, :kmax]), -np.amin(psi_2[:ymax, :kmax]))
    #     max = np.amax(psi_2[:ymax, :kmax])
    #     a = ax2.contourf(Y,Z,psi_2[:ymax, :kmax], levels=np.linspace(-max,max,1e2), extend = 'both')
    #     a.cmap.set_under('yellow')
    #     a.cmap.set_over('cyan')
    #     plt.colorbar(a, ax=ax2)
    #     ax3 = axes[0,2]
    #     ax3.set_title('psi_4')
    #     max = np.maximum(np.amax(psi_4[:ymax, :kmax]), -np.amin(psi_4[:ymax, :kmax]))
    #     a = ax3.contourf(Y,Z,psi_4[:ymax, :kmax], levels=np.linspace(-max,max,1e2))
    #     plt.colorbar(a, ax=ax3)
    #
    #     ax = axes[1,1]
    #     ax.set_title('w')
    #     var = w[ic,:ymax,:kmax]
    #     max = np.maximum(np.amax(var), -np.amin(var))
    #     a = ax.contourf(Y,Z,var[:ymax, :kmax], levels=np.linspace(-max,max,1e2), cmap=cm_bwr)
    #     speed_yz = np.sqrt(v[ic,:ymax,:kmax]*v[ic,:ymax,:kmax] + var*var)
    #     lw = 5 * speed_yz[:, :] / speed_yz[:, :].max()
    #     ax.streamplot(y,z,v[ic,:ymax,:kmax].T,var.T, color='k', linewidth=lw[:,:].T)
    #     ax3.streamplot(y,z,v[ic,:ymax,:kmax].T,var.T, color='k', linewidth=lw[:,:].T)
    #     plt.colorbar(a, ax=ax)
    #     ax = axes[1,2]
    #     ax.set_title('w')
    #     max = np.maximum(np.amax(psi_4[:ymax, :kmax]), -np.amin(psi_4[:ymax, :kmax]))
    #     a = ax.contourf(Y, Z, psi_4[:ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax)
    #     var = w[ic,:,:]
    #     max = np.maximum(np.amax(var[:ymax, :kmax]), -np.amin(var[:ymax, :kmax]))
    #     a = ax.contour(Y,Z,var[:ymax, :kmax], levels=np.linspace(-max,max,1e1+1), cmap=cm_bwr)
    #     # plt.colorbar(a, ax=ax)
    #
    #     ax = axes[1,0]
    #     ax.set_title('s')
    #     var = s[ic,:,:]
    #     max = np.amax(var[:ymax, :kmax])
    #     min = np.amin(var[:ymax, :kmax])
    #     a = ax.contourf(Y,Z,var[:ymax, :kmax], levels=np.linspace(min, max, 1e2))
    #     plt.colorbar(a, ax=ax)
    #
    #     fig.savefig(os.path.join(path_out_figs, 'psi_t'+str(t0)+'s.png'))
    #     plt.close(fig)
    #
    #     fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    #     plt.suptitle('t='+str(t0)+'s')
    #     ax1 = axes[0]
    #     ax1.set_title('v')
    #     max = np.maximum(np.amax(v[ic,20:20 + ymax, :kmax]), -np.amin(v[ic,20:20 + ymax, :kmax]))
    #     a = ax1.contourf(Y, Z, v[ic, 20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax1 = axes[1]
    #     ax1.set_title('w')
    #     max = np.maximum(np.amax(w[ic,20:20 + ymax, :kmax]), -np.amin(w[ic,20:20 + ymax, :kmax]))
    #     a = ax1.contourf(Y, Z, w[ic, 20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax2 = axes[2]
    #     ax2.set_title('vort_yz')
    #     print vort_yz.shape
    #     max = np.maximum(np.amax(vort_yz[20:20 + ymax, :kmax]), -np.amin(vort_yz[20:20 + ymax, :kmax]))
    #     a = ax2.contourf(Y, Z, vort_yz[20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax2)
    #     ax2 = axes[3]
    #     ax2.set_title('psi_4')
    #     max = np.maximum(np.amax(psi_4[20:20 + ymax, :kmax]), -np.amin(psi_4[20:20 + ymax, :kmax]))
    #     a = ax2.contourf(Y, Z, psi_4[20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax2)
    #     fig.savefig(os.path.join(path_out_figs, 'psi_2_t' + str(t0) + 's.png'))
    #     plt.close(fig)
    #
    #
    #     ym = np.int(3*ny_/3)
    #     ym = ny
    #     YM_, ZM_ = np.meshgrid(y_arr[:ym],z_arr[:kmax],indexing='ij')
    #     # compare streamfunction and streamlines
    #     fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))
    #     # # plt.imshow(psi[:,:40].T, origin='lower')
    #     a = ax1.contourf(YM_,ZM_,psi_1[:ym,:kmax])
    #     plt.colorbar(a, ax=ax1)
    #     # y_arr = np.arange(0,ny)
    #     # z_arr = np.arange(0,ny)
    #     ax1.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T, density=1.5, color='k')
    #     # print '.....', y_arr.shape, z_arr.shape, ym, kmax, YM_.shape, ZM_.shape, v[ic,:ym,:kmax].shape, w[ic,:ym,:kmax].shape
    #     # ax1.streamplot(YM_,ZM_,v[ic,:ym,:kmax], w[ic,:ym,:kmax], density=1.5, color='k')
    #     b = ax2.contourf(y_arr[:ym], z_arr[:kmax], vort_yz[:ym,:kmax].T)
    #     ax2.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T)
    #     plt.colorbar(b, ax=ax2)
    #     ax1.set_title('stream function (x(i-1)..x(i+1)')
    #     ax2.set_title('vorticity')
    #     plt.suptitle('t='+str(t0)+'s')
    #     plt.savefig(os.path.join(path_out_figs, 'psi_1_t'+str(t0)+'s.png'))
    #     plt.close()
    #
    #     print('')

    return


# ----------------------------------------------------------------------

def streamfunction_circular_flow():
    print '--- compute circular stream ---'

    x = np.arange(-10,10,0.1)
    y = np.arange(-10,10,0.1)


    X, Y = np.meshgrid(x, y, indexing='ij')  # if indexing not defined, output will be 'xy'-indexed

    U = Y / (X**2 + Y**2)
    V = - X / (X**2 + Y**2)

    umax = np.amax(U)
    vmax = np.amax(V)
    fig, axes = plt.subplots(2,4, figsize=(20,10))
    ax1 = axes[0,0]
    a = ax1.contourf(x,y,U.T, levels=np.linspace(-umax,umax), cmap=cm_bwr)
    ax1.set_title('contourf(x,y,U.T)')
    ax = axes[1,0]
    ax.contourf(X,Y,U, levels=np.linspace(-umax,umax), cmap=cm_bwr)
    ax.set_title('contoruf(X,Y,U)')
    plt.colorbar(a, ax=ax1)
    plt.colorbar(a, ax = ax)
    ax2 = axes[0,1]
    b = ax2.contourf(x,y,V.T, levels=np.linspace(-10, 10), cmap=cm_bwr)
    ax2.set_title('contourf(x,y,V.T)')
    ax = axes[1,1]
    ax.contourf(X,Y,V, levels=np.linspace(-10,10), cmap=cm_bwr)
    ax.set_title('contourf(X,Y,V)')
    plt.colorbar(b, ax = ax2)
    plt.colorbar(b, ax = ax)

    ax3 = axes[0,2]
    x_ = np.arange(-10,10,0.8)
    y_ = np.arange(-10,10,0.8)
    X_, Y_ = np.meshgrid(x_, y_)
    U_ = Y_ / (X_ ** 2 + Y_ ** 2)
    V_ = - X_ / (X_ ** 2 + Y_ ** 2)
    # q = ax3.quiver(x_, y_, X_, Y_)
    # ax3.quiverkey(q, X=0.3, Y=1.1, U=10,
    #              label='Quiver key, length = 10', labelpos='E')
    q = ax3.quiver(x_, y_, U_, V_)
    ax3.quiverkey(q, X=0.3, Y=1.1, U=3,
                label='Quiver key, length = 3', labelpos='E')
    ax3.set_title('quiver(x,y,U,V)')
    ax = axes[1,2]
    q = ax.quiver(X_, Y_, U_, V_)
    ax.set_title('quiver(X,Y,U,V)')
    # q = ax3.quiver(x, y, U.T, V.T)
    ax.quiverkey(q, X=0.3, Y=1.1, U=3,
                 label='Quiver key, length = 3', labelpos='E')

    # psi = np.log(X**2 + Y**2)
    # a = ax4.contourf(x, y, psi)
    psi_ = np.log(X_**2 + Y_**2)
    ax4 = axes[0,3]
    a = ax4.contourf(x_, y_, psi_)
    plt.colorbar(a, ax=ax4)
    ax4.streamplot(x_, y_, U_, V_)
    # CS = ax4.contour(x, y, psi)
    # ax4.clabel(CS)
    ax4.set_title('streamplot(x,y,U,V)')
    ax = axes[1,3]
    ax.contourf(x_, y_, psi_)
    plt.colorbar(a, ax=ax)
    ax.streamplot(X_,Y_,U_,V_)
    ax.set_title('streamplot(X,Y,U,V)')
    fig.savefig('./circ_flow.png')
    plt.close()

    return


def compute_streamfunction_simps5(v_, w_, rho0, i0, j0):
    # int dy w(y,0)
    yrange = np.arange(w_.shape[0])
    # psi_A_ = integrate.simps(w_[:, 0], yrange[:])
    rhow = np.zeros(shape=w_.shape)
    for k in range(nz):
        rhow[:,k] = rho0[k]*w_[:,k]
    psi_A = integrate.simps(w_[:,:], yrange[:], axis=0)
    print 'shapes: ', psi_A.shape
    # print 'diff: ', (psi_A[0]-psi_A_)
    krange = np.arange(v_.shape[1])
    psi_B = integrate.simps(v_[0,:], krange[:], axis=0)

    return psi_A


def compute_streamfunction_simps4(v_, w_, rho0, t0):
    # Psi(y,z) = \int_y0^y w(y',z)\rho_0(z) dy' + a(z)
    # a(z) = -\int_z0^z u(y0,z) dz'
    if v_.shape != w_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [ly,lz] = v_.shape
    psi_A = np.zeros((ly,lz), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((ly,lz), dtype=np.double)     # psi_B = int dy u(x,y)
    # psi_simps = np.zeros((ly,lz))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for j in range(ly):
        for k in range(lz):
            # Simpson Rule
            j0 = np.int(np.round(j/2))
            k0 = np.int(np.round(k/2))
            # (1)  \rho_0(z)*\int_y0^y w(y',z) dy'
            psi_A[j,k] = rho0[k] * (j * dy) / 6 * (w_[0, k] + 4 * w_[j0, k] + w_[j, k])
            # (2) a(z) =  -\int_z0^z u(y0,z) dz'
            psi_B[j,k] = -  j * dy / 6 * (rho0[0]*v_[0, 0] + 4 * rho0[k0]*v_[0, k0] + rho0[k]*v_[0, k])
    psi_simps = psi_A + psi_B

    return psi_simps


def compute_streamfunction_simps3(u_, v_, t0):
    # Psi(x,y) = \int_y0^y w(x,y') dy' + a(x)
    # a(x) = -\int_x0^x v(x',y0) dx'
    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [lx, ly] = u_.shape
    psi_A = np.zeros((lx, ly), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((lx, ly), dtype=np.double)     # psi_B = int dy u(x,y)
    psi_simps = np.zeros((lx, ly))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for i in range(lx):
        for j in range(ly):
            # Simpson Rule
            i0 = np.int(np.round(i/2))
            j0 = np.int(np.round(j/2))
            # (1)  \int_y0^y u(x,y') dy'
            psi_B[i, j] = j * dy / 6 * (u_[i, 0] + 4 * u_[i, j0] + u_[i, j])
            # (2) a(x) = -\int_x0^y v(x',y0)
            psi_A[i, j] = - i * dx / 6 * (v_[0, 0] + 4 * v_[i0, 0] + v_[i, 0])
    psi_simps = psi_A + psi_B

    return psi_simps


def compute_streamfunction_simps2(u_, v_, t0):
    # F(x) = int_c^x f(y) dy, with F(c) = 0
    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [lx, ly] = u_.shape
    psi_A = np.zeros((lx, ly), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((lx, ly), dtype=np.double)     # psi_B = int dy u(x,y)
    psi_simps = np.zeros((lx, ly))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for i in range(lx):
        for j in range(ly):
            # Simpson Rule
            i0 = np.int(np.round(i/2))
            j0 = np.int(np.round(j/2))
            psi_B[i, j] = j * dy / 6 * (u_[i, 0] + 4 * u_[i, j0] + u_[i, j])
            psi_A[i, j] = i * dx / 6 * (v_[0, j] + 4 * v_[i0, j] + v_[i, j])
    psi_simps = psi_A - psi_B

    return psi_simps


def compute_streamfunction_simps(u_, v_, t0):
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
    y_arr = dy*np.arange(0,ny_)
    z_arr = dz*np.arange(0,nz)
    YM_, ZM_ = np.meshgrid(y_arr[:ym], z_arr[:km])
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
    fig.savefig(os.path.join(path_out_figs, 'test_u_t'+str(t0)+'s.png'))
    plt.close()

    return psi_simps


def integrand(x, a, b):
    return a*x**2 + b

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

    global path, path_fields, path_out_figs, path_out_fields
    path = args.path
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_vorticity')
    path_out_fields = os.path.join(path, 'fields_vorticity')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    if not os.path.exists(path_out_fields):
        os.mkdir(path_out_fields)

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]
    print('nx, ny, nz', nx, ny, nz)

    ''' determine file range '''
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
        kmax = 60

    print('times', timerange)
    print('kmax ', kmax)

    return timerange, kmax, nml
