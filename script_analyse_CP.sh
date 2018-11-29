#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

do_loop=0
dTh=$1

echo "dTh:" $dTh

path="/cond1/meyerbe/ColdPools/single_3D_noise/"
casename="ColdPoolDry_single_3D"

if [ $dTh -eq 4 ]
then 
  z_params=( 1730, 870, 430 )
  r_params=( 430, 870, 1730 )
elif [ $dTh -eq 3 ]
then
  #z_params=( 2000 500 1000 2000 )
  #r_params=( 2000 2000 1000 500 )
  z_params=( 2000 1000 )
  r_params=( 2000 1000 )
elif [ $dTh -eq 2 ]
then
  z_params=( 2450 1225 815 )
  r_params=( 815 1225 2450 )
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
echo "z-parameters:" $z_params 
echo "r-parameters:" $r_params

echo "#geometry parameters:" $n_geom

echo "MIN MAX ALL"
python compute_minmax_all.py $casename $path --zparams ${z_params[*]} --rparams ${r_params[*]}



#count_geom=0
#while [ $count_geom -lt $n_geom ]
#do
#  zstar=${z_params[$count_geom]}
#  rstar=${r_params[$count_geom]}
#  echo "parameters:" $zstar $rstar
#  
#  id="dTh"$dTh"_z"$zstar"_r"$rstar
#  echo $id
#  
#  fullpath=$path$id
#  echo $fullpath
  
#  # PLOT STREAMLINES 
#  #python plot_streamlines_singleCP.py ColdPoolDry_single_3D $fullpath
#  # MIN MAX VALUES (w, s)
#  python compute_minmax.py $casename $fullpath --tmin 100 --tmax 3600

#  echo " " 
#  ((count_geom++))
#done
#
#fi

#for i in ${th_params[@]}; do
#  echo $i
#done

echo "finished bash script"

