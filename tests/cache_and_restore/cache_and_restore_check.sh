rm plot_log*

make sim_test

time ./sim_test abc_covid.json --simulate --serial 0

CKSUM=`cksum plot_log_*`
NUM_MATCH=`cksum plot_log_* | cut -d' ' -f1,2 | uniq -c | wc -l`

if [[ $NUM_MATCH -eq 1 ]]
then
    echo "CACHE AND RESTORE TEST PASSED"
    echo "ALL CHECKSUMS EQUIVALENT"
    echo "$CKSUM"
else
    echo "CACHE AND RESTORE TEST FAILED"
    echo "$CKSUM"
fi
