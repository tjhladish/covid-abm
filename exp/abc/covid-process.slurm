#!/bin/bash
#SBATCH --job-name=covid_abc
#SBATCH --output=./auto_output/covid_%A_%a.out
#SBATCH --error=./auto_output/covid_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --chdir=/home/tjhladish/work/covid-abm/exp/abc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=0:10:00
#SBATCH --array=0
#SBATCH --partition=hpg2-compute

module load gcc/7.3.0 gsl

./abc_sql abc_covid.json --process 2>> fl_v5a.0.log
