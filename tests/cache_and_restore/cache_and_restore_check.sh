#!/bin/bash

rm plot_log*

SERIAL=$1

if [[ $# -eq 0 ]]
then
    echo "ERROR: pass sim serial as only argument"
    exit 1
fi

time ./sim_test abc_covid.json --simulate --serial $SERIAL

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
