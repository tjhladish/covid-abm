#!/bin/bash
#SBATCH --job-name=covid_abc
#SBATCH --output=./auto_output/covid_%A_%a.out
#SBATCH --error=./auto_output/covid_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --chdir=/home/tjhladish/work/covid-abm/exp/state-vac
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=2:00:00
#SBATCH --array=0-999
##SBATCH --partition=default

module load gcc gsl

for i in `seq 1 5`;
do
    ./sim_test abc-covid-active-v4.0.json --simulate
done
