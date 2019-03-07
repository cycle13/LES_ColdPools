import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib.pyplot as plt


def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--k0")
    parser.add_argument("--var_name")
    args = parser.parse_args()

    times, nml = set_input_parameters(args)
    dt_stats = nml['stats_io']['frequency']
    print('stats output timestep: ', dt_stats)
    print ''


    # read in stats file
    nt = len(times)
    variable = np.ndarray((nt, nz))


    var_list = ['viscosity_mean']
    var_list = ['u_max', 'v_max', 'w_max']
    var_list = ['s_min', 'theta_mean']
    var_list = ['dynamic_pressure_mean']
    if args.var_name:
        var_list = [args.var_name]

    for var_name in var_list:
        print '-------'+var_name + '------------'
        restart = 1
        count = 0
        itmax = 0
        while restart:
            if count == 0:
                file_name = 'Stats.'+case_name + '.nc'
            else:
                file_name = 'Stats.' + case_name + '.Restart_' + str(count-1) + '.nc'
            print file_name
            try:
                stats_file_ = nc.Dataset(os.path.join(path, 'stats', file_name))
            except:
                try:
                    stats_file_ = nc.Dataset(os.path.join(path, file_name))
                except:
                    restart = 0
            count += 1

            stats_time = stats_file_.groups['timeseries'].variables['t'][:]
            ts_grp = stats_file_.groups['timeseries'].variables
            prf_grp = stats_file_.groups['profiles'].variables


            if restart:
                print 'time: ', stats_time, count
                print variable.shape

                stats_var = prf_grp[var_name][:,:]

                for it, t0 in enumerate(stats_time):
                    print '--t', t0, it
                    if t0 <= times[-1]:
                        it_tot = np.int(np.double(t0)/dt_stats) - 1
                        print it_tot
                        variable[it_tot, :] = stats_var[it,:]
                # itmax = stats_time[-1]/dt_stats
                # print 'itmax', itmax

        stats_file_.close()



        # plotting
        if var_name == 'diffusivity_mean':
            max = 0.8
        elif var_name == 'w_max':
            max = 3.5
        elif var_name == 'v_max':
            max = 10.
        elif var_name == 's_min':
            max = 6871
        elif var_name == 'theta_min':
            max = 300.2

        if var_name == 's_min':
            min = 6860
        elif var_name == 'theta_min':
            min = 296.9
        else:
            min = 0

        fig_name = var_name + '.png'
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey='all', figsize=(11, 5))
        # lvls = np.linspace(np.amin(viscosity_mean), np.amax(viscosity_mean), 1e3)
        lvls = np.linspace(min, max, 1e3)
        cf = ax1.contourf(times, np.arange(0, kmax)*dx[2], variable[:, :kmax].T, levels=lvls)
        plt.colorbar(cf, ax=ax1)
        ax1.set_title(var_name)
        lvls = np.linspace(np.amin(variable), np.amax(variable), 1e3)
        cf = ax2.contourf(times, np.arange(0, kmax)*dx[2], variable[:, :kmax].T, levels=lvls)
        plt.colorbar(cf, ax=ax2)
        ax1.set_title(var_name)
        ax1.set_xlabel('time')
        ax2.set_xlabel('time')
        plt.suptitle('dz='+str(dx[2])+'m')
        plt.tight_layout
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return


# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_fields, path_out_data, path_out_figs
    path = args.path
    path_fields = os.path.join(path, 'fields')

    path_out_data = os.path.join(path, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path, 'figs_stats')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    print path_out_data
    print path_out_figs
    print path_fields
    print ''

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']

    global kmax
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nz

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    ''' time range '''
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()

    print('')
    print('times', times)
    print('')
    print('kmax ', kmax, 'nx ', nx)
    print('')

    return times, nml
# _______________________________________________________
# _______________________________________________________
def theta_s(s):
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# _______________________________________________________

if __name__ == '__main__':
    main()
