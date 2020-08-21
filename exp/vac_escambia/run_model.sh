for i in `seq 1 360`
do
./abc_sql abc_covid.json --simulate 2> /dev/null
done
