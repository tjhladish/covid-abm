reps=200
seedstart=$((1+$1*$reps))
seedend=$(($reps+$1*$reps))

echo $seedstart ' to ' $seedend

for s in `seq 0 3`
do
    for r in `seq $seedstart $seedend`
    do
        ./sim_test $s $r 2>> mmod_$1
    done
done
