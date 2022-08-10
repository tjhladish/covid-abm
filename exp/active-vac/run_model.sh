for i in `seq 1 120`
do
    ./sim_test abc-flvtms-vac_comp.json --simulate 2> /dev/null > /dev/null
done
