jobid=`sbatch covid-simulate.slurm | grep -P -o '\d+'`

for i in `seq 0 10`
do
jobid=`sbatch --depend=afterany:$jobid covid-process.slurm | grep -P -o '\d+'`
jobid=`sbatch --depend=afterany:$jobid covid-simulate.slurm | grep -P -o '\d+'`
done

jobid=`sbatch --depend=afterany:$jobid covid-process.slurm | grep -P -o '\d+'`
