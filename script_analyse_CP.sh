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


path="/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2/"
casename="ColdPoolDry_single_3D"

# to run over different dTh, use input dTh=0
if [ $dTh -eq 0 ]
then
  do_r1km=1
  #dTh_params=( 2 3 4 )
else
  do_r1km=0
  #dTh_params=( $dTh )
fi
#echo "dTh_params: " ${dTh_params[@]}


# set geometry parameters
if [ $dTh -eq 1 ]
then
  z_params=( 3465 1730 1155 )
  r_params=( 1155 1730 3465 )
elif [ $dTh -eq 2 ]
then
  #z_params=( 2450 1225 815 ) #run1
  #r_params=( 815 1225 2450 ) #run1 
  z_params=( 500 900 1600 1900 2500 ) #run2
  r_params=( 1900 1300 900 800 600 )  #run2
elif [ $dTh -eq 3 ]
then
  #z_params=( 2000 250 500 670 1000 1500 2000 4000 ) #run1
  #r_params=( 2000 4000 2000 1500 1000 670 500 250 ) #run1
  z_params=( 500 1000 1600 2000 2500 ) #run2
  r_params=( 1500 1000 700 600 500 )   #run2
  z_params=( 500 1000 2000 ) #run3
  r_params=( 1500 1000 600 ) #run3
  z_params=( 1000 )
  r_params=( 1000 )

elif [ $dTh -eq 4 ]
then
  #z_params=( 430 870 1730 ) #run1
  #r_params=( 1730 870 430 ) #run1
  z_params=( 500 900 1600 2000 2500 ) #run2
  r_params=( 1300 900 600 500 400 )   #run2
  #z_params=( 2500 )
  #r_params=( 400 )
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
echo " "



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

  echo "make smaller files"
  #python convert_fields_smaller_k.py $casename $fullpath --kmax 80
  #python convert_fields_smaller_k.py $casename $fullpath --vert True

  #echo "MIN MAX VALUES (w, s)"
  #python compute_minmax.py $casename $fullpath --tmin 100 --tmax 3600

  #echo "CROSSSECTIONS"
  #python plot_crosssections.py $casename $fullpath --tmin 100 --tmax 800


  #echo "CP HEIGHT"
  #python plot_CP_height.py $casename $fullpath --tmin $tmin --tmax $tmax
  
  #echo "PLOT STREAMLINES"
  #python plot_streamlines_singleCP.py ColdPoolDry_single_3D $fullpath --tmin $tmin --tmax $tmax
  
  #echo "CP RIM"
  ## for each simulation compute CP rim (a) Bettina, (b) Marielle
  ##     >> r(phi,t), U_r(phi,t), r_av(phi,t), U_r,av(phi,t)
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --kmin 0 --kmax 5
  #python define_cp_rim_nbi_v2.py $casename $fullpath --tmin $tmin --tmax $tmax --perc 80
  ## >> compare radius, radial velocity r(phi,t), U_r(phi,t); average radius; r_av(t); average rim velocity U_r,av(t)

  ##echo "ENERGY"
  ##python compute_energy_domain.py $casename $fullpath --tmin $tmin --tmax $tmax

  #echo "VORTICITY"
  #python vorticity_streamfct_compute.py $casename $fullpath --tmin $tmin --tmax $tmax

  echo " "
  ((count_geom++))
done

#echo "CP height all"
#python plot_CP_height_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "plot CP RIM all"
#python plot_CP_rim_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

echo "plot CP RIM from tracers"
z_params_r1km=( 900 1000 900 )
r_params_r1km=( 1300 1000 900 )
python plot_tracer_analysis_all.py ColdPoolDry_single_3D $path $dTh --zparams_r1km ${z_params_r1km[*]} --rparams_r1km ${r_params_r1km[*]} --k0 0 --tmin 100 --tmax 3600

#echo "ENERGY all"
#python compute_energy_domain_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax

#echo "VORTICITY all"
#python vorticity_streamfct_plotting_all.py $casename $path $dTh --zparams ${z_params[*]} --rparams ${r_params[*]} --tmin $tmin --tmax $tmax





# -------------------------------------------
do_r1km=0

if [ $do_r1km -eq 1 ]
then
# r=1km
dTh_params=( 2 3 4 )
z_params=( 900 1000 900 )
r_params=( 1300 1000 900 )
echo "dTh-parameters:" ${dTh_params[@]}
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}


#echo "MIN MAX r1km"
#python compute_minmax_r1km.py $casename $path ${dTh_params[*]} ${z_params[*]} ${r_params[*]} --tmin $tmin --tmax $tmax
#echo " "


fi #end r=1km

echo "finished bash script"

