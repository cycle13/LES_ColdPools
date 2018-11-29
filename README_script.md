



# Cold Pools
**Plotting Tools**
> all python files can be called directly or through bash-script `script_analyse_CP.sh`

***(a) plot initial configuration:***
`plot_configuration.py` [casename CASENAME] [path_root PATH] [dTh DELTA_TH]
[--zparams Z_RANGE] [--rparams R_RANGE] [--tmin TMIN] [--tmax TMAX]

> OUTPUT: figures in 'path_root/figs_config/'

***(b) plot min/max in domain and xy-crosssections:*** `compute_minmax_all.py`

***(c) plot CP height:*** `plot_CP_height.py`
[casename CASENAME] [path_root PATH][--s_crit SCRIT] 
[--tmin TMIN] [--tmax TMAX][--kmin KMIN] [--kmax KMAX]

computes CP height in terms of entropy threshold (given by s_crit)
> OUTPUT: plots with CP height, 
max vertical velocity and height of maximum vertical velocity






