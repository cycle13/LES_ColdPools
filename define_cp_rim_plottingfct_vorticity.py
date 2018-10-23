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
