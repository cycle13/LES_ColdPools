#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)


dTh=$1

path="/cond1/meyerbe/ColdPools/single_3D_noise/"

if [ $dTh -eq 3 ]
then
  z_params=( 2000 500 1000 2000 )
  r_params=( 2000 2000 1000 500 )
elif [ $dTh -eq 2 ]
then
  z_params=( 2450 1225 815 )
  r_params=( 815 1225 2450 )
elif [ $dTh -eq 1 ]
then 
  z_params=( 3465 1730 1155 )
  r_params=( 1155 1730 3465 )
fi


n_geom=${#z_params[@]}
n_therm=${#th_params[@]}
n_tot=$(( $n_geom*$n_therm ))

echo "dTh:" $dTh
echo "#geometry parameters:" $n_geom



count_geom=0
while [ $count_geom -lt $n_geom ]
do
  # 'echo $z_params' will output only the first element.
  zstar=${z_params[$count_geom]}
  rstar=${r_params[$count_geom]}
  echo "parameters:" $zstar $rstar
  
  id="dTh"$dTh"_z"$zstar"_r"$rstar
  echo $id
  
  fullpath=$path$id
  echo $fullpath
  
  python plot_streamlines_singleCP.py ColdPoolDry_single_3D $fullpath
  

  # use the sleep command to add delay fora  specified amoutn of time
  # s for seconds (default); m for minutes; h for hours; d for days
  #sleep 10

  echo " " 
  ((count_geom++))
done

#for i in ${th_params[@]}; do
#  echo $i
#done

echo "finished bash script"

