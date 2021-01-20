#!/bin/bash
#SBATCH --job-name=covid_abc
#SBATCH --output=./auto_output/covid_%A_%a.out
#SBATCH --error=./auto_output/covid_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --chdir=/home/tjhladish/work/covid-abm/exp/abc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --time=12:00:00
#SBATCH --array=0-199
#SBATCH --partition=hpg2-compute

module load gcc/7.3.0 gsl

for i in `seq 1 60`;
do
    ./abc_sql abc_covid.json --simulate -n 1
done
