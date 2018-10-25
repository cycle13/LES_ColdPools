import numpy as np
import matplotlib.pyplot as plt
import os

# global cm_bwr, cm_grey, cm_vir
# cm_bwr = plt.cm.get_cmap('bwr')
# cm_vir = plt.cm.get_cmap('viridis')
# cm_grey = plt.cm.get_cmap('gist_gray_r')

def set_colorbars(cm_bwr_, cm_vir_, cm_grey_):

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = cm_bwr_
    cm_grey = cm_grey_
    cm_vir = cm_vir_
    return


def plot_vorticity(v, w, vort, t0, path):
    kmax = 50

    plt.figure(figsize=(15,20))
    plt.subplot(311)
    plt.imshow(v[:,:kmax].T, origin='lower')
    plt.title('v')
    plt.colorbar(shrink=0.5)
    plt.subplot(312)
    plt.imshow(w[:,:kmax].T, origin='lower')
    plt.title('w')
    plt.colorbar(shrink=0.5)
    plt.subplot(313)
    plt.imshow(vort[:,:kmax].T, origin='lower')
    plt.title('vort_yz')
    plt.colorbar(shrink=0.5)
    plt.savefig(os.path.join(path, 'vort_yz_t'+str(t0)+'.png'))
    plt.close()

    print(os.path.join(path, 'vort_yz_t'+str(t0)+'.png'))

    return

def plot_histogram(data, data_name, perc_range, t0, path):
    kmax = 50

    n_bins = 1e2
    min = np.amin(data)
    max = np.amax(data)
    bins = np.linspace(min, max, n_bins)

    perc = np.zeros(len(perc_range))
    for i,p in enumerate(perc_range):
        perc[i] = np.percentile(data[:,:,:kmax], p)

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))
    plt.subplot(211)
    n, bins, patches = plt.hist(data[:, :, :kmax].flatten(), bins)
    max_hist = np.amax(n)
    for i, p in enumerate(perc_range):
        plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=1, label=str(p) + 'th percentile')
    plt.ylabel('n')
    plt.legend(loc='best')  # , bbox_to_anchor=(0, 0))

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    n_bins = 50
    bins_pos = np.linspace(0, max, n_bins)
    logbins = np.logspace(np.log10(bins_pos[1]), np.log10(bins_pos[-1]), len(bins_pos))
    plt.subplot(212)
    plt.hist(data[:, :, :kmax].flatten(), bins=logbins)
    plt.xscale('log')
    max_hist = 1e3
    for i, p in enumerate(perc_range):
        plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=3, label=str(p) + 'th percentile')
    plt.legend(loc='right')  # , bbox_to_anchor=(0, 0))

    plt.xlabel(data_name)
    plt.suptitle('Histogram of ' + data_name + '   (t=' + str(t0) + r',  k$\in$[0,' + str(kmax) + '])')
    plt.savefig(os.path.join(path, 'hist_'+data_name+'_kmax' + str(kmax) + '_t' + str(t0) + '.png'))
    plt.close()


    return