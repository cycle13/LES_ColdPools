import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc



def main():
    ampl = 10
    res = '1km'
    path_data = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/333-600_T0_300K_ampl_' \
                +str(ampl)+'_' + res + '/IRT_out'

    path_data1_d1 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/1-288_T0_300K_ampl_10_1km/IRT_out'
    path_data1_d2 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/333-600_T0_300K_ampl_10_1km/IRT_out'
    path_data1_d3 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/621-888_T0_300K_ampl_10_1km/IRT_out'
    path_data1_d4 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/909-1176_T0_300K_ampl_10_1km/IRT_out'
    root1 = 'T0_300K_ampl_10_1km'
    path_data2 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/333-600_T0_300K_ampl_10_1km/IRT_out'
    path_data3 = '/nbi/ac/conv2/haerter/cold_pool_tracking/cp_tracking_Feb2019/333-600_T0_300K_ampl_10_1km/IRT_out'
    filename = 'irt_tracks_output_pure_r_int.txt'
    path_figs = '/nbi/ac/cond1/meyerbe/checkerboard/Iorg'

    f = open(os.path.join(path_data, filename), "r")
    contents =f.read()
    f.close()

    f1 = np.loadtxt(os.path.join(path_data1_d1, filename), skiprows=0)
    f2 = np.loadtxt(os.path.join(path_data1_d2, filename), skiprows=0)
    f3 = np.loadtxt(os.path.join(path_data1_d3, filename), skiprows=0)
    f4 = np.loadtxt(os.path.join(path_data1_d4, filename), skiprows=0)

    # plot COM
    figname = root1 + 'COM.png'
    fig, axes = plt.subplots(1,4, figsize=(20, 6))
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]
    ax1.plot(f1[:,8], f1[:,9], 'ko', markersize=1)
    ax2.plot(f2[:,8], f2[:,9], 'ko', markersize=1)
    ax3.plot(f3[:,8], f3[:,9], 'ko', markersize=1)
    ax4.plot(f4[:,8], f4[:,9], 'ko', markersize=1)
    fig.suptitle('COM')
    for ax in axes.flat:
        ax.set_aspect('equal')
    fig.savefig(os.path.join(path_figs, figname))
    plt.close(fig)


    # I_org
    f = f4
    del f1, f2, f3, f4
    # loop through timesteps
    # for each timestep and each cell compute NN-distance

    # get timesteps: (column 2)
    n_time = np.int(f[-1,1])
    print('number of timestesp: ' + str(n_time))
    n_time = 2

    # get range of cell-IDs
    max_ID = np.int(np.max(f[:,0]))
    print('max ID: ' + str(max_ID))
    # max_ID = 100


    min_dist = 999999.9*np.ones((max_ID+1, n_time), dtype=np.double)
    i_start_cells = np.zeros(max_ID+2, dtype=np.int)
    l = 0
    for ID in range(1,max_ID+1):
        l = np.min(np.where(f[:,0] == ID))
        i_start_cells[ID] = l
        print(ID, 'l', l)
    print('')

    # max_d = 15. # search-radius: 15km = 15pts
    max_d = 400. # search-radius: 15km = 15pts
    for ID in range(1, max_ID + 1):
        l = i_start_cells[ID]
        it = 0
        t = np.int(f[l+it,1])
        print('--- ID', ID, l, t)
        x = f[l,8]
        y = f[l,9]
        for ID_ in range(1, max_ID+1):
            if ID_ != ID:
                l_1 = i_start_cells[ID_]
                l_2 = i_start_cells[ID_+1]
                t_ = np.int(f[l_1,1])
                # print('ID_', f[l_1,0], l_1, t_, l_2, np.int(f[l_2-1,1]))
                if np.int(f[l_1,1])<=t and np.int(f[l_2-1,1]) >= t:
                    while f[l_1,0] == ID_ and (t_ < t):
                        ## print('l_, t_', l_1, t_, f[l_1,0])
                        l_1 += 1
                        t_ = np.int(f[l_1,1])
                    # print('t_', t_, t)
                    while (t_ == t) and (f[l_1,0]==ID_):
                        # print('testing', ID_, t_, np.int(f[l_1,1]), np.abs(f[l_1,8] - x), np.abs(f[l_1,9] - y))
                        if np.abs(f[l_1,8] - x) <= max_d and np.abs(f[l_1,9] - y) <= max_d:
                            d = np.sqrt((f[l_1,8]-x)**2 + (f[l_1,9]-y)**2)
                            # print('dist: ',  d, min_dist[ID,it])
                            min_dist[ID,it] = np.minimum(min_dist[ID,it], d)
                        l_1 += 1
                        t_ = f[l_1, 1]
    print(min_dist[:,0])
    dump_file('min_dist.nc', path_figs, n_time, max_ID, min_dist)



    figname = root1 + '_NN.png'
    data = [i for i in min_dist[:,0] if i < 999999.9]
    n_bins = 1e2
    min = 0
    max = np.amax(data)
    print('min', min, 'max', max)
    # perc = np.zeros(len(perc_range))
    # for i, p in enumerate(perc_range):
    #     perc[i] = np.percentile(data[:, :, :kmax], p)

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))
    bins = np.linspace(min, max, 1e2)
    n, bins, patches = axes[0].hist(min_dist[:, 0].flatten(), bins)
    # max_hist = np.amax(n)
    # print('n', np.amax(n), n.shape)
    bins = np.linspace(min, max, 5e1)
    n, bins, patches = axes[1].hist(data, bins, normed=False, facecolor='red', alpha=0.4, label='w')
    # for i, p in enumerate(perc_range):
    #     plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=1, label=str(p) + 'th percentile')

    # # histogram on log scale.
    # # Use non-equal bin sizes, such that they look equal on log scale.
    # n_bins = 50
    # bins_pos = np.linspace(0, max, n_bins)
    # print('bins', bins)
    # # logbins = np.logspace(np.log10(bins_pos[1]), np.log10(bins_pos[-1]), len(bins_pos))
    # # axes[0].hist(data[:].flatten(), bins=logbins)
    # axes[1].hist(data[:], bins=bins_pos)
    # # plt.xscale('log')
    # # max_hist = 4e4
    # # for i, p in enumerate(perc_range):
    # #     plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=3, label=str(p) + 'th percentile')
    # # plt.legend(loc='right')  # , bbox_to_anchor=(0, 0))

    for ax in axes.flatten():
        ax.set_xlabel('NND [km]')
        ax.set_ylabel('n')
    # plt.ylabel('log_10(n)')
    # plt.suptitle('Histogram of w' + '   (t=' + str(t0) + r',  k$\in$[0,' + str(kmax) + '])')
    plt.savefig(os.path.join(path_figs, figname))
    plt.close()




    return

#_______________________________

# def dump_file(fname, dTh_range, zstar_range, rstar_range,
#               PE_ref, scaling, PE_ref_range, PE_numerical,
#               params_dict, count, path_out):

def dump_file(fname, path_out, n_time, n_cells, min_dist):
    rootgrp = nc.Dataset(os.path.join(path_out, fname), 'w', format='NETCDF4')
    rootgrp.createDimension('nt', n_time)
    rootgrp.createDimension('n_cells', n_cells+1)
    var = rootgrp.createVariable('NND', 'f8', ('n_cells', 'nt'))
    print('dumping', min_dist.shape, n_time, n_cells)
    var[:,:] = min_dist[:,:]
    # var = rootgrp.createVariable('rstar', 'f8', ('n_rstar'))
    # var[:] = rstar_range[:]
    # var = rootgrp.createVariable('dTh', 'f8', ('n_dTh'))
    # var[:] = dTh_range[:]
    # var = rootgrp.createVariable('scaling', 'f8', ('n_scal'))
    # var[:] = scaling[:]
    # var = rootgrp.createVariable('PE_ref', 'f8', )
    # var[:] = PE_ref
    # var = rootgrp.createVariable('PE_ref_scaled', 'f8', ('n_scal'))
    # var[:] = PE_ref_range[:]
    # var = rootgrp.createVariable('PE_numerical', 'f8', ('n_dTh', 'n_rstar'))
    # var[:,:] = PE_numerical[:,:]
    #
    # for i, s in enumerate(scaling):
    #     var = rootgrp.createVariable('parameters_s'+str(s), 'f8', ('n', 'params'))
    #     n = params_dict[str(s)].shape[0]-1
    #     print('s', s, params_dict[str(s)].shape)
    #     if n > 0:
    #         var[:n,:] = params_dict[str(s)][1:,1:]
    rootgrp.close()
    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
