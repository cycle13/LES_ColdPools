



# Cold Pools
**Plotting Tools**
> all python files can be called directly or through bash-script `script_analyse_CP.sh`

***(a) plot initial configuration:***

***`plot_configuration.py`*** [casename CASENAME] [path_root PATH] [dTh DELTA_TH]
[--zparams Z_RANGE] [--rparams R_RANGE] [--tmin TMIN] [--tmax TMAX]

> OUTPUT: figures in 'path_root/figs_config/' 
> - contourf xy-plane of entropy
> - envelopes of initial anomaly (cos2())


***(b) plot min/max in domain and xy-crosssections:*** 

***`compute_minmax_all.py`***

***(c) plot CP height:*** 

***`plot_CP_height.py`***
[casename CASENAME] [path_root PATH][--s_crit SCRIT] 
[--tmin TMIN] [--tmax TMAX][--kmin KMIN] [--kmax KMAX]

computes CP height in terms of entropy threshold (given by s_crit)
> OUTPUT: 
> - plots with CP height in *id/figs_CP_height/*, 
max vertical velocity and height of maximum vertical velocity
> - nc-files with max(CP_height[t]), max(w)[i,j], height(max(w))[i,j]

***`plot_CP_height_all.py`***: reads in output fields from `plot_CP_height.py` and plots it 
for all CPs to compare


***(d) plot streamlines:*** 

***`plot_streamlines.py`, `plot_streamlines_singleCP.py`***
[casename CASENAME] [path PATH]
[--tmin TMIN] [--tmax TMAX][--kmin KMIN] [--kmax KMAX]

plots crosssection with contourfigure of variable (`cont_var`) overlaid by 
streamlines of 2D-velocity field (projection onto crosssection plane); streamlines
computed by`matplotlib.pylab.streamplot`, thickness (and color) weighted by 
absolute value of velocity 

`plot_streamplot_xy`,`plot_streamplot_xz`, `plot_streamplot_yz`
> OUTPUT: figures in 'path/streamlines/' 

 





