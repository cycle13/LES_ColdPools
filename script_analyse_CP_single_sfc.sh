#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
# read -p "dTh: " dTh; 
dTh=3;
read -p "tmin: " tmin; 
read -p "tmax: " tmax;
read -p "kmax: " kmax; 

# if no input for tmin, tmax
if [[ $tmin = "" ]]
  then
  tmin=100
fi
if [[ $tmax = "" ]]
  then 
  tmax=3600
fi
echo "tmin=" $tmin ", tmax=" $tmax

if [[ $kmax = "" ]]
  then 
  kmax=120
fi
echo "kmax=" $kmax
echo" "

path="/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_on/single_3D_noise/"
casename="ColdPoolDry_single_3D"

z_params=( 1000 ) 
r_params=( 1000 )
dx_params=( 25 50 100 ) 

n_geom=${#r_params[@]}
#n_therm=${#th_params[@]}
n_res=${#dx_params[@]}
#n_tot=$(( $n_geom*$n_therm ))
echo "dTh:" $dTh
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}
echo "res-parameters:" ${dx_params[@]}
echo "#resolution parameters:" $n_res

echo " "
#echo "TEST INITIALIZATION / CONFIGURATION"
#python plot_configuration.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]}
echo " "


count_res=0
while [ $count_res -lt $n_res ]
do
  zstar=${z_params[0]}
  rstar=${r_params[0]}
  dx=${dx_params[$count_res]}
  echo "parameters:" $zstar $rstar
  
  id="dTh"$dTh"_z"$zstar"_r"$rstar"_dx"$dx
  echo $id
  
  fullpath=$path$id
  echo $fullpath
  echo " "

  #echo "make smaller files"
  #python convert_fields_smaller_k.py $casename $fullpath --vert True --hor True --kmax 80 --tmin $tmin --tmax $tmax
  #python convert_fields_smaller_k.py $casename $fullpath --vert True --tmin $tmin --tmax $tmax

  #echo "CP HEIGHT"
  #python CP_height_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

  echo "average angular"
  python average_angular.py $casename $fullpath --tmin $tmin --tmax $tmax --kmax $kmax 

  echo " "
  ((count_res++))
done


#echo "plot CP RIM from tracers"
#python plot_tracer_analysis_all_PEscaling.py $casename $path $dTh --k0 0 --tmin $tmin --tmax $tmax

#echo "MIN MAX ALL"
#python compute_minmax_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax
#echo " "

#echo "CP height all"
#python plot_CP_height_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM all"
#python plot_CP_rim_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM from tracers"
#z_params_r1km=( 900 1000 900 )
#r_params_r1km=( 1300 1000 900 )
#python plot_tracer_analysis_all.py ColdPoolDry_single_3D $path $dTh --zparams_r1km ${z_params_r1km[*]} --rparams_r1km ${r_params_r1km[*]} --k0 0 --tmin 100 --tmax 3600

#echo "ENERGY all"
#python compute_energy_domain_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "VORTICITY all"
#python vorticity_streamfct_plotting_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax





# -------------------------------------------

echo "finished bash script"

