#!bin/bash/
#python tracers_plot_2D_singleCP_dx50m.py
#    --tmin 0 --tmax 3600 --timerange 100 900 2100 --imin 50 --shift_x -0.5 --shift_t 0

shift_t=( 0 1 )
shift_x=-0.5

t0_range=( 0 100 )
t1_range=( 900 )
t2_range=( 1800 2100 2400 2700 )
timerange=( $t0 $t1 $t2 )

#python tracers_plot_2D_singleCP_dx50m.py --tmin 0 --tmax 3600 --timerange 100 200 300 --imin 50 --shift_x -0.5 --shift_t 0

echo " "
#echo ${shift_t[@]}
for st in ${shift_t[@]}
do
    for t0 in ${t0_range[@]}
    do
        for t1 in ${t1_range[@]}
        do
            for t2 in ${t2_range[@]}
            do
                echo "shift_t " $st ", (t0, t1, t2) " $t0 $t1 $t2
                if [ $t2 -lt 2200 ]
                then
                    imin=50
                else
                    imin=0
                fi
                python tracers_plot_2D_singleCP_dx50m.py --tmin 0 --tmax 3600 --timerange $t0 $t1 $t2 --imin $imin --shift_x -0.5 --shift_t $st
                echo " "
            done
        done

    done
done