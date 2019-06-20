#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
read -p "dTh: " dTh; 
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
  tmax=100
fi
echo "tmin=" $tmin ", tmax=" $tmax

if [[ $kmax = "" ]]
  then
  kmax=120
fi
echo "kmax=" $kmax
echo" "


path="/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/"
casename="ColdPoolDry_single_3D"

# to run over different dTh, use input dTh=0
if [ $dTh -eq 0 ]
then
  do_r1km=1
else
  do_r1km=0
fi
#echo "dTh_params: " ${dTh_params[@]}


# set geometry parameters
if [ $dTh -eq 1 ]
then
  z_params=( 3465 1730 1155 )
  r_params=( 1155 1730 3465 )
elif [ $dTh -eq 2 ]
then
  z_params=( 2500 1900 1600 900 500 )
  r_params=( 600 800 900 1300 1900 )
elif [ $dTh -eq 3 ]
then
  z_params=( 2500 2000 1600 1000 500 )
  r_params=( 500 600 700 1000 1500 )
  z_params=( 2500 2000 1000 500 )
  r_params=( 500 600 1000 1500 )
elif [ $dTh -eq 4 ]
then 
  z_params=( 2500 2000 900 500 )
  r_params=( 400 500 900 1300 )
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
  #python convert_fields_smaller_k.py $casename $fullpath --vert True --tmin $tmin --tmax $tmax
  #python convert_fields_smaller_k.py $casename $fullpath --hor True --tmin $tmin --tmax $tmax

  #echo "MIN MAX VALUES (w, s)"
  #python compute_minmax.py $casename $fullpath --tmin 100 --tmax 3600

  #echo "CROSSSECTIONS"
  #python plot_crosssections.py $casename $fullpath --tmin 100 --tmax 800

  #echo "compute CP HEIGHT"
  #python CP_height_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

#  echo "ANGULAR AVERAGE"
#  python average_angular.py $casename $fullpath --kmax 20 --tmin $tmin --tmax $tmax

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
  #python vorticity_streamfct_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

  echo " "
  ((count_geom++))
done


#echo "MIN MAX ALL"
#python compute_minmax_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax
#echo " "

echo "CP height all"
python CP_height_plot_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM all"
#python plot_CP_rim_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM from tracers"
#z_params_r1km=( 900 1000 900 )
#r_params_r1km=( 1300 1000 900 )
#python tracer_analysis_all.py $casename $path $dTh --zparams_r1km ${z_params_r1km[*]} --rparams_r1km ${r_params_r1km[*]} --k0 0 --tmin 100 --tmax 3600

#echo "ENERGY all"
#python compute_energy_domain_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "VORTICITY all"
#python vorticity_streamfct_plotting_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax


# -------------------------------------------

echo "finished bash script"

