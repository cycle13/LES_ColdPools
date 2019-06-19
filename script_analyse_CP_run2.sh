#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
read -p "dTh: " dTh; 
read -p "tmin: " tmin; 
read -p "tmax: " tmax; 

# if no input for tmin, tmax
if [[ $tmin = "" ]]
  then
  tmin=100
fi
if [[ $tmax = "" ]]
  then 
  tmax=100
fi

echo "tmin=" $tmin ", tmax=" $tmax


path="/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/"
casename="ColdPoolDry_single_3D"

if [ $dTh -eq 4 ]
then 
  z_params=( 2500 2000 900 500 )
  r_params=( 400 500 900 1300 )
elif [ $dTh -eq 3 ]
then
  z_params=( 2500 2000 1600 1000 500 )
  r_params=( 500 600 700 1000 1500 )
  #z_params=( 4000 )
  #r_params=( 250 )
elif [ $dTh -eq 2 ]
then
  z_params=( 2500 1900 1600 900 500 )
  r_params=( 600 800 900 1300 1900 )
  #z_params=( 815 )
  #r_params=( 2450 )
elif [ $dTh -eq 1 ]
then 
  z_params=( 3465 1730 1155 )
  r_params=( 1155 1730 3465 )
elif [ $dTh -eq 10 ]
then 
  z_params=( 2000 )
  r_params=( 2000 )
fi


n_geom=${#z_params[@]}
n_therm=${#th_params[@]}
n_tot=$(( $n_geom*$n_therm ))
echo "dTh:" $dTh
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}

echo "#geometry parameters:" $n_geom




echo " "

#echo "TEST INITIALIZATION / CONFIGURATION"
#python plot_configuration.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]}
#echo " "

#echo "MIN MAX ALL"
#python compute_minmax_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax
#echo " "


count_geom=0
while [ $count_geom -lt $n_geom ]
do
  zstar=${z_params[$count_geom]}
  rstar=${r_params[$count_geom]}
  echo "parameters:" $zstar $rstar
  
  id="dTh"$dTh"_z"$zstar"_r"$rstar
  echo $id
  
  fullpath=$path$id
  echo $fullpath
  echo " "

  #echo "make smaller files"
  #python convert_fields_smaller_k.py $casename $fullpath --kmax 80

  #echo "MIN MAX VALUES (w, s)"
  #python compute_minmax.py $casename $fullpath --tmin 100 --tmax 3600

  #echo "ANGULAR AVERAGE"
  #python average_angular.py $casename $fullpath --kmax 10 --tmin $tmin --tmax $tmax

  echo "compute CP HEIGHT"
  python CP_height_compute.py $casename $fullpath --tmin $tmin --tmax $tmax
  
  #echo "PLOT STREAMLINES"
  #python plot_streamlines_singleCP.py ColdPoolDry_single_3D $fullpath --tmin $tmin --tmax $tmax
  
  #echo "CP RIM"
  ## for each simulation compute CP rim (a) Bettina, (b) Marielle
  ##     >> r(phi,t), U_r(phi,t), r_av(phi,t), U_r,av(phi,t)
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --kmin 0 --kmax 5
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --perc 80
  ## >> compare radius, radial velocity r(phi,t), U_r(phi,t); average radius; r_av(t); average rim velocity U_r,av(t)

  #echo "ENERGY"
  #python compute_energy_domain.py $casename $fullpath --tmin $tmin --tmax $tmax

  #echo "VORTICITY"
  #python vorticity_streamfct_compute.py --casename $casename --path $fullpath --tmin $tmin --tmax $tmax

  echo " "
  ((count_geom++))
done

#echo "CP height all"
#python plot_CP_height_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM all"
#python plot_CP_rim_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin 100 --tmax 300

# echo "plot CP RIM from tracers"
#plot_tracer_analysis_all.py ColdPoolDry_single_3D /cond1/meyerbe/ColdPools/3D_sfc_fluxes_on/single_3D_noise/ 2 --k0 0 --tmin 100 --tmax 3600

#echo "ENERGY all"
#python compute_energy_domain_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "VORTICITY all"
#python vorticity_streamfct_plotting_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

echo "finished bash script"

