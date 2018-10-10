



# Cold Pools

**(1) Definition of Cold Pool rim**
>> Define rim of cold pool by a threshold (usually 95th percentile) of the vertical velocity

(1a) `define_cp_rim.py`: 
Find outer rim of mask based on number of neighbours
INPUT: [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--k0 K0]

OUTPUT:
  - figures in ``PATH/figs_cp_rim``

(1b) `define_cp_rim_v2.py`: Find inner and outer rim of mask based on number of neighbours and filling interior of mask
INPUT: [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--kmin KMIN] [--kmax KMAX]

OUTPUT:
- figures in ``PATH/figs_cp_rim``  
- 3D fields ``PATH/fields_cp_rim/rimmask_perc**th_t**.nc``
    - mask[x,y,k] (x=0..nx_, y=0..ny_, k=kmin..kmax)
    - rim_inner[x,y,k], rim_outer[x,y,k]
    - profiles:
        - krange[k], zrange[k] (k=kmin..kmax)
        - k_dumped[k] \in {0,1} (0=not dumped, 1=dumped) 
 
 
(1c) `define_cp_rim_v3.py`: like v2 but with changes to boundaries to improve performance     

define_cp_rim_plottingfct.py


Details:
(1a) `define_cp_rim.py`: Find outer rim of mask based on number of neighbours
INPUT: [--casename CASENAME] [--path PATH] [--tmin TMIN] [--tmax TMAX] [--k0 K0]
1. read in w-field, shift field (roll) and define partial domain where to look for cold pool
2. mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
3. Define rim of cold pool as the outline of the mask; based on number of neighbours


**HOW TO USE:**

`CloudClosure.py`:

_input:_ reads in nc-files with 3D field output, computes (1) the bivariate Gaussian from maximum likelihood estimation and (2) the empirical PDF from Kernel Density Estimation (KDE)

_output:_ PDF parameters for Gaussian PDF




**Kernel Density Estimation**
The principle behind nearest _neighbor methods_ is to find a predefined number of training samples closest in distance to the new point, and predict the label from these.

The implemented Kernel Density algorithm uses the Ball Tree or KD Tree for efficient queries

- *Kernel Density Estimation*: `from sklearn.neighbors import KernelDensity`
`class sklearn.neighbors.KernelDensity(bandwidth=1.0, algorithm='auto', kernel='gaussian', metric='euclidean', atol=0, rtol=0, breadth_first=True, leaf_size=40, metric_params=None)`
- parameters:
    - kernel = 'gaussian' (‘gaussian’|’tophat’|’epanechnikov’|’exponential’|’linear’|’cosine’; default is 'gaussian')
    - bandwith = number (controlls "smoothness" of estimation)
- methods:
    - `KernelDensity.fit(data)`
    - `KernelDensity.fit(data).score_samples()`:  Evaluate the density model on the data.
    - `KernelDensity.fit(data).score(data)`:          total log probability
    - `KernelDensity.fit(data).sample(data)`: generate random variables from model




