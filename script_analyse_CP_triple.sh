#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
#read -p "dTh: " dTh;
#dTh=3
dTh=5
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

path="/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/"
casename="ColdPoolDry_triple_3D"

# set geometry parameters
if [ $dTh -eq 3 ]
then
  z_params=( 2000 )
  r_params=( 2000 )
  d_params=( 10 15 20 )
elif [ $dTh -eq 5 ]
then 
  z_params=( 1000 )
  r_params=( 1100 ) 
  d_params=( 10 15 20 )
fi


n_geom=${#d_params[@]}
n_therm=${#th_params[@]}
n_tot=$(( $n_geom*$n_therm ))
echo "dTh:" $dTh
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}
echo "#geometry parameters:" $n_geom
echo " "
echo "tmin=" $tmin ", tmax=" $tmax
echo " "

#echo "TEST INITIALIZATION / CONFIGURATION"
#python plot_configuration.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]}
#echo " "


count_geom=0
while [ $count_geom -lt $n_geom ]
do
  zstar=${z_params[0]}
  rstar=${r_params[0]}
  d=${d_params[$count_geom]}
  echo "parameters:" $zstar $rstar $d
  
  id="dTh"$dTh"_z"$zstar"_r"$rstar"_d"$d"km"
  echo $id
  
  fullpath=$path$id
  echo $fullpath
  echo " "

#  echo "make smaller files"
#  python convert_fields_smaller_k.py $casename $fullpath --vert True --tmin $tmin --tmax $tmax --kmax $kmax
#  python convert_fields_smaller_k.py $casename $fullpath --hor True --tmin $tmin --tmax $tmax --kmax 3
#
#  echo "MIN MAX VALUES (w, s)"
#  python compute_minmax_tripleCP.py $casename $fullpath --tmin $tmin --tmax $tmax
#  echo " "

#  echo "CROSSSECTIONS"
#  python plot_crosssections.py $casename $fullpath --tmin $tmin --tmax $tmax
#  echo " "
#
#  echo "CP HEIGHT"
#  python CP_height_compute.py $casename $fullpath --tmin $tmin --tmax $tmax
#  echo " "

  #echo "PLOT STREAMLINES"
  #python plot_streamlines_triple.py $casename $fullpath --tmin $tmin --tmax $tmax --hor True

  #echo "CP RIM"
  ## for each simulation compute CP rim (a) Bettina, (b) Marielle
  ##     >> r(phi,t), U_r(phi,t), r_av(phi,t), U_r,av(phi,t)
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --kmin 0 --kmax 5
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --perc 80
  ## >> compare radius, radial velocity r(phi,t), U_r(phi,t); average radius; r_av(t); average rim velocity U_r,av(t)

  ##echo "ENERGY"
  ##python compute_energy_domain.py $casename $fullpath --tmin $tmin --tmax $tmax

#  echo "VORTICITY"
#  python vorticity_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

  echo " "
  ((count_geom++))
done

echo "CP height all"
python plot_CP_height_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

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







echo "finished bash script"

