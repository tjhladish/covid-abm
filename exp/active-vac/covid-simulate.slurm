#!/bin/bash
#SBATCH --job-name=covid_abc
#SBATCH --output=./auto_output/covid_%A_%a.out
#SBATCH --error=./auto_output/covid_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --chdir=/home/tjhladish/work/covid-abm/exp/active-vac
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-199
##SBATCH --partition=default

module load gcc gsl

sleep `echo ${SLURM_ARRAY_TASK_ID}%200 | bc`
for i in `seq 1 750`;
do
    #./sim_test3.0 abc-par_update-fl3.0.json --simulate --doFlorida
    #./sim_test abc-ppb_fit-v4-w_season.json --simulate
    #./sim_test abc-covid-active-v6.1.json --simulate
    ./sim_test abc-covid-active-v7.json --simulate
done
