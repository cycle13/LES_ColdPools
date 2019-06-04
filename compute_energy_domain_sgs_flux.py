import numpy as np
import netCDF4 as nc
import os

''' copied algorithms from PyCLEs, momentum_diffusion.h '''

def compute_velocity_gradient(d_advecting, d_advected, nx, ny, nk, dx, t0, path_fields):
    # print('compute velocity gradient', d_advecting, d_advected)
    vgrad = np.zeros((nx, ny, nk), dtype=np.double)
    dxi = 1./dx[d_advecting]
    if d_advected == 0:
        vel_name = 'u'
    elif d_advected == 1:
        vel_name = 'v'
    elif d_advected == 2:
        vel_name = 'w'
    root = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'), 'r')
    vel = root.groups['fields'].variables[vel_name][:,:,:nk]
    root.close()

    if d_advecting == 0:
        for i in range(nx-1):
            vgrad[i, :, :] = (vel[i+1,:,:] - vel[i,:,:]) * dxi;
    elif d_advecting == 1:
        for j in range(ny-1):
            vgrad[:, j, :] = (vel[:,j+1,:] - vel[:,j,:]) * dxi;
    elif d_advecting == 2:
        for k in range(nk-1):
                vgrad[:, :, k] = (vel[:,:,k+1] - vel[:,:,k]) * dxi;
    del vel
    return vgrad




def compute_strain_rate(d_advecting, d_advected, nx, ny, nk, dx, t0, path_fields):
    # print('compute strain rate', d_advecting, d_advected)
    # strain_rate = np.zeros((nx, ny, nk))
    vgrad_1 = compute_velocity_gradient(d_advecting, d_advected, nx, ny, nk, dx, t0, path_fields)
    if d_advected != d_advecting:
        vgrad_2 = compute_velocity_gradient(d_advected, d_advecting, nx, ny, nk, dx, t0, path_fields)
        strain_rate = 0.5 * (vgrad_1 + vgrad_2)
        del vgrad_2
    else:
        strain_rate = np.array(vgrad_1, copy=True)
    del vgrad_1

    return strain_rate




def compute_diffusive_flux(viscosity_mean, rho0, d_advecting, d_advected, nx, ny, nk, dx, t0, path_fields):
    # print('compute diffusive flux', d_advecting, d_advected)
    flux = np.zeros((nx, ny, nk), dtype = np.double)
    strain_rate = compute_strain_rate(d_advecting, d_advected, nx, ny, nk, dx, t0, path_fields)

    # Compute flux if not flux of $\tau_{3,3}$ or $\tau{3,2}$, or $\tau{3,1}$, $\tau{2,3}$, or $\tau{1,3}$
    if d_advecting != 2 and d_advected != 2:
        # Viscosities are located at the flux location so no interpolation necessary
        for k in range(nk-1):
            flux[:,:,k] = -2.0 * strain_rate[:,:,k] * viscosity_mean[k+1] * rho0[k]
    # Compute flux if $\tau_{3, 3}$
    elif d_advecting == 2 and d_advected == 2:
        # Viscosity again at flux location so no interpolation required
        for k in range(nk-1):
            flux[:,:,k] = -2.0 * strain_rate[:,:,k] * viscosity_mean[k+1] * rho0[k+1]
    # Compute flux if $\tau_{3, 1}$, $\tau_{3, 2}$, $\tau_{1, 3}$, $\tau_{2, 3}
    else:
        # Viscosity requires interpolation, but since viscosity_mean, only interpolation in k-direction
        for k in range(nk-1):
            visc_interp = 0.5 * (viscosity_mean[k] + viscosity_mean[k + 1])
            flux[:,:,k] = -2.0 * strain_rate[:,:,k] * visc_interp * rho0[k]

    return flux

