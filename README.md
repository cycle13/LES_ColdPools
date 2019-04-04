



# Cold Pools

#### (1) Definition of Cold Pool rim:
*Define rim of cold pool by a threshold (usually 95th percentile) of the vertical velocity*

**(1a) Find outer rim of mask based on number of neighbours**

***`define_cp_rim.py`*** [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--k0 K0] 
> OUTPUT: figures in ``PATH/figs_cp_rim``

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


Details:
(1a) `define_cp_rim.py`: Find outer rim of mask based on number of neighbours
INPUT: [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--k0 K0]
1. read in w-field, shift field (roll) and define partial domain where to look for cold pool
2. mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
3. Define rim of cold pool as the outline of the mask; based on number of neighbours


#### (2) Compute Energy Budget

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




#### (3) Vorticity computation CP rim
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
### ***averaging etc.***

#### (1) Fields of reduced Dimension:
***`convert_fields_smaller_k.py`*** [path PATH] [--kmin KMIN] [--kmax KMAX] [--k0 K0]

*make new field of reduced number of levels (k=kmin..kmax) or single level (k=k0) and merged time*

***modules:***

- ``convert_file_for_varlist.py``: one single output file of n levels (k=kmin..kmax) for all variables in 
given var_list and for all times

- ``convert_file_for_varlist_vertsection.py``: one single output file with yz_crosssection for all variables in 
given var_list and for all times

- ``convert_file_for_singlevariable_onelevel.py:`` returns one file with all timesteps at one level (k=k0) for 
given variable 
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