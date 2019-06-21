#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
dTh=5;
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


path="/nbi/ac/conv3/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/"
casename="ColdPoolDry_single_3D"

z_params=( 1000 )
r_params=( 500 1100 1600 2300 )

n_geom=${#r_params[@]}
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
  zstar=${z_params[0]}
  rstar=${r_params[$count_geom]}
  echo "parameters:" $zstar $rstar
  
  id="dTh"$dTh"_z"$zstar"_r"$rstar
  echo $id
  
  fullpath=$path$id
  echo $fullpath
  echo " "

  #echo "make smaller files"
  #python convert_fields_smaller_k.py $casename $fullpath --vert True --tmin $tmin --tmax $tmax --kmax $kmax
  #python convert_fields_smaller_k.py $casename $fullpath --hor True --tmin $tmin --tmax $tmax  --kmax $kmax

  #echo "MIN MAX VALUES (w, s)"
  #python compute_minmax.py $casename $fullpath --tmin 100 --tmax 3600

  #echo "CROSSSECTIONS"
  #python plot_crosssections.py $casename $fullpath --tmin 100 --tmax 800

  #echo "compute CP HEIGHT"
  #python CP_height_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

  echo "compute CP VOLUME"
  python CP_volume_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

#  echo "ANGULAR AVERAGE"
#  python average_angular.py $casename $fullpath --tmin $tmin --tmax $tmax --kmax 20

  #echo "PLOT STREAMLINES"
  #python plot_streamlines_singleCP.py ColdPoolDry_single_3D $fullpath --tmin $tmin --tmax $tmax


  echo " "
  ((count_geom++))
done


#echo "MIN MAX ALL"
#python compute_minmax_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax
#echo " "

echo "CP height all"
python CP_height_volume_plot_all.py $casename $path $dTh --rmax_plot 9000 --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM all"
#python plot_CP_rim_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM from tracers"
#z_params_r1km=( 900 1000 900 )
#r_params_r1km=( 1300 1000 900 )
#echo $path
#python tracer_analysis_all_PEscaling.py $casename $path $dTh --k0 0 --tmin $tmin --tmax $tmax


#echo "ENERGY all"
#python compute_energy_domain_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "VORTICITY all"
#python vorticity_streamfct_plotting_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax


# -------------------------------------------

echo "finished bash script"

