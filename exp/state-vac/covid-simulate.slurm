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
#SBATCH --mem-per-cpu=20000
#SBATCH --time=10:00:00
#SBATCH --array=0-399
##SBATCH --partition=default

module load gcc gsl

for i in `seq 1 1`;
do
    #./sim_test abc-covid-states-v2.0.json --simulate
    ./sim_test-no_db abc-covid-states-v2.0.json --simulate --serial ${SLURM_ARRAY_TASK_ID}
done
