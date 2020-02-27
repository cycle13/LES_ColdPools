import matplotlib.pyplot as plt

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['text.usetex'] = 'true'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1.
plt.rcParams['ytick.major.width'] = 1.
plt.rcParams['ytick.minor.width'] = 1.
# plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['lines.linewidth'] = 2
# plt.rcParams['grid.linewidth'] = 20
plt.rcParams['pdf.fonttype'] = 42         # Output Type 3 (Type3) or Type 42 (TrueType)
plt.rcParams['font.sans-serif'] = 'Helvetica'


cm = plt.cm.get_cmap('rainbow')
# cm_bwr = plt.cm.get_cmap('bwr')
cm_bwr = plt.cm.get_cmap('seismic')
cm_bwr_r = plt.cm.get_cmap('seismic_r')
cm_bw = plt.cm.get_cmap('bone_r')
cm_bw_r = plt.cm.get_cmap('bone')
cm_gray = plt.cm.get_cmap('gray')
cm_contourfigs = cm_bwr
# cm_twin = plt.cm.get_cmap('twilight')

colorlist_all = ['darkred', 'maroon', 'r', 'tomato', 'indianred', 'orange', 'gold',
                 'limegreen', 'forestgreen', 'g', 'darkgreen', 'seagreen', 'lightseagreen',
                 'darkcyan',
                 'mediumblue', 'darkblue', 'midnightblue', 'navy']
colorlist = ['maroon', 'indianred', 'orange',
                     'limegreen', 'darkgreen', 'darkcyan', 'lightseagreen',
                     'mediumblue', 'navy']
colorlist5 = ['maroon', 'indianred', 'orange', 'darkcyan', 'navy']
    # colorlist5 = ['orange', 'indianred', 'maroon', 'navy', 'lightseagreen']
colorlist4 = ['maroon', 'orange', 'darkcyan', 'navy']
colorlist4_blue = ['orange', 'green', 'dodgerblue', 'mediumblue']
colorlist3 = ['indianred', 'darkcyan', 'navy']
# colorlist2 = ['maroon', 'navy']
colorlist2 = ['maroon', 'tomato']

colorlist2_bw = [cm_bw(0.5), cm_bw(1.)]

# ----------------------------------------------------------------------

'''time windows for collisions (dTh=5K, z*=1000m)'''
# rstar = 1100m
d_range_r1100 = [10,    12,     15]
t_ini_r1100 =   [600,   600,    600]
t_2CP_r1100 =   [1100,  1500,   2400]
# t_3CP_r1100 =   [1500,  2200,   3300]
t_3CP_r1100 =   [1600,  2200,   3300]
t_final_r1100 = [2500,  3100,   4200]

# rstar = 2000m
d_range_r2000 = [10,    15,     20]
t_ini_r2000 =   [400,   400,    400]
t_2CP_r2000 =   [600,   800,    2200]
t_3CP_r2000 =   [800,   1100,   3200]
t_final_r2000 = [1800,  2500,   3500]

d_range = {}
d_range['1100'] = d_range_r1100
d_range['2000'] = d_range_r2000
t_ini = {}
t_ini['1100'] = t_ini_r1100
t_ini['2000'] = t_ini_r2000
t_2CP = {}
t_2CP['1100'] = t_2CP_r1100
t_2CP['2000'] = t_2CP_r2000
t_3CP = {}
t_3CP['1100'] = t_3CP_r1100
t_3CP['2000'] = t_3CP_r2000
t_final = {}
t_final['1100'] = t_final_r1100
t_final['2000'] = t_final_r2000