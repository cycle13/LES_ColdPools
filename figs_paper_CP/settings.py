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
                 'limegreen', 'forestgreen', 'g', 'darkgreen', 'seagreen', 'lightseagreen', 'darkcyan',
                 'mediumblue', 'darkblue', 'midnightblue', 'navy']
colorlist = ['maroon', 'indianred', 'orange',
                     'limegreen', 'darkgreen', 'darkcyan', 'lightseagreen',
                     'mediumblue', 'navy']
colorlist5 = ['maroon', 'indianred', 'orange', 'darkcyan', 'navy']
    # colorlist5 = ['orange', 'indianred', 'maroon', 'navy', 'lightseagreen']
colorlist4 = ['indianred', 'orange', 'darkcyan', 'navy']

colorlist2 = ['maroon', 'navy']
colorlist2 = [cm_bw(0.5), cm_bw(1.)]
