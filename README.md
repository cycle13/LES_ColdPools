



# Cold Pools

#### (1) Definition of Cold Pool rim:
*Define rim of cold pool by a threshold (usually 95th percentile) of the vertical velocity*

**(1a) Find outer rim of mask based on number of neighbours**

***`define_cp_rim.py`*** [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--k0 K0] 
> OUTPUT: figures in ``PATH/figs_cp_rim``

Details:
1. read in w-field, shift field (roll) and define partial domain where to look for cold pool
2. mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
3. Define rim of cold pool as the outline of the mask; based on number of neighbours


**(1b) Find inner and outer rim of mask based on number of neighbours and filling interior of mask**

***`define_cp_rim_v2.py`*** [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--kmin KMIN] [--kmax KMAX]
 
> OUTPUT:
> - figures in ``PATH/figs_cp_rim``  
> - 3D fields ``PATH/fields_cp_rim/rimmask_perc**th_t**.nc``
>    - mask[x,y,k] (x=0..nx_, y=0..ny_, k=kmin..kmax)
>    - rim_inner[x,y,k], rim_outer[x,y,k]
>    - profiles:
>        - krange[k], zrange[k] (k=kmin..kmax)
>        - k_dumped[k] \in {0,1} (0=not dumped, 1=dumped) 
 
 
**(1c) like v2 but with changes to boundaries to improve performance**

`define_cp_rim_v3.py`: 

`define_cp_rim_plottingfct.py`




#### (2) Definition of Cold Pool rim width:
*Computes CP rim width based on dipole of vertical updraft w: distance between positions of maximum and 
minimum vertical velocity at a given level (width=r(max(w))-r(min(w))) 

***`rim_width.py`***
[--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--kmin KMIN][--kmax KMAX]

> OUTPUT:
> - figures in ``PATH/figs_radial_average``  

#### (3) Computation of Cold Pool Height:
*Computes the CP height based on threshold on entropy anomaly for each column, starting from the top*
***`plot_CP_height.py`*** [casename CASENAME] [path PATH] [--tmin TMIN] [--tmax TMAX] [--kmax KMAX] [--s_crit S_CRIT]
> OUTPUT: stats-file and figures in ``PATH/figs_CP_height``
>
> stats-file: 
> - 2D-fields (for each time step) for CP-height, max(w) and height(max(w)) for each column
> - timeseries with domain maximum values of CP-height, max(w) and height(max(w))


 #### (4) Vorticity computation CP rim
***`conceptual_model/vorticity_rim.py`*** [casename CASENAME][path PATH]

*computes and saves 4D field (k=0,..kmax) for radial velocity and azimuthal average of all variables*

 > OUTPUT: 
 > - 4D field of radial velocity v_r(t, x, y, z) with z=0,..kmax*dz; filename: `path/fields_v_rad/v_rad.nc`
 > - stats-file with data of form var(t,r,z) in `path/data_analysis/stats_radial_averaged.nc`
 
 ***`average_angular_plot_compare.py`*** [casename CASENAME]
 [path_root PATH_ROOT][path1 PATH1][path2 PATH][path_out PATH_OUT][--tmin TMIN][--tmax TMAX][--kmax KMAX]
 (Note: path1, path2 and path_out are relative paths that are concatenated with path_root)
 *


### _______________________________________
### ***Energy***
#### (1) Compute Energy Budget

> (a) compute kinetic energy (KE)
>
> KE ~ v^2 = (u^2 + v^2 + w^2)
>
> (b) compute potential energy (PE) by integrating over anomaly
>
> PE = \int dz g * (th_anomaly(z) - th_env(z)) * z

(2a) ```compute_pot_energy_v1.py```
> compute KE and PE in area that is defined by masked on w 
(same idea as for definition of CP rim, but freshly computed)

(2b) ```compute_pot_energy_v2.py```
> compute KE and PE in area that is defined by mask, 
computed in definition of CP rim (read in file)

***module ``compute_KE.py``***

***module ``compute_PE.py``*** 






### _______________________________________
### ***averaging etc.***

#### (1) Fields of reduced Dimension:
***`convert_fields_smaller_k.py`*** [casename CASENAME] [path PATH] [--tmin TMIN] [--tmax TMAX]
[--kmin KMIN] [--kmax KMAX] [--k0 K0][--k0 K0] [--vert BOOL] [--hor BOOL]

*make new field of reduced number of levels (k=kmin..kmax) or single level (k=k0) and merged time*

***modules:***

- ``convert_file_for_varlist.py``: one single output file of n levels (k=kmin..kmax) for all variables in 
given var_list and for all times

- ``convert_file_for_varlist_vertsection.py``: one single output file with yz_crosssection for all variables in 
given var_list and for all times

- ``convert_file_for_varlist_horsection.py``: one single output file with xy_crosssection (at level ``k0``) 
for all variables in given var_list and for all times

- ``convert_file_for_singlevariable_onelevel.py:`` returns one file with all timesteps at one level (k=k0) for 
single given variable 
(used for Olga's tracers) 



#### (2) Azimuthal Average 
***`average_angular.py`*** [casename CASENAME][path PATH][--tmin TMIN][--tmax TMAX][--kmax KMAX]

*computes and saves 4D field (k=0,..kmax) for radial velocity and azimuthal average of all variables*

 > OUTPUT: 
 > - 4D field of radial velocity v_r(t, x, y, z) with z=0,..kmax*dz; filename: `path/fields_v_rad/v_rad.nc`
 > - stats-file with data of form var(t,r,z) in `path/data_analysis/stats_radial_averaged.nc`
 
 ***`average_angular_plot_compare.py`*** [casename CASENAME]
 [path_root PATH_ROOT][path1 PATH1][path2 PATH][path_out PATH_OUT][--tmin TMIN][--tmax TMAX][--kmax KMAX]
 (Note: path1, path2 and path_out are relative paths that are concatenated with path_root)
 *
 
### _______________________________________
### ***Initialization***
#### (1) Test PyCLES Initialization
***`Initialization_singleCP.py`***

***`Initialization_tripleCP.py`***

#### (2) Initialization of same Potential Energy:

***`Initialization_PE.py`***

> OUTPUT: 
> - nc-files for each dTh (2, 3, 4K) with z*- and r*-values that correspond to the same reference potential energy
(PE[dTh=3K, r*=1km, z*=1km])  


### _______________________________________
### ***Thermodynamics***
#### (1) Thermodynamic functions from PyCLES:
***`thermodynamic_functions.py`***