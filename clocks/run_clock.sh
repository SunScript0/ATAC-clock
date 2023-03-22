#!/bin/sh
#SBATCH --partition=standard
#SBATCH --ntasks=11
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2G
#SBATCH -o tpm_all_samples2.log
#SBATCH -t 10:00:00
#SBATCH --exclude=bhc0062,bhc0035

path="/gpfs/fs2/scratch/fmorandi/ChromAcc-clock/clocks/parallel/2023-03-19_10-06_tpm_all_samples"
srun="srun --exclusive --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"

for i in {0..10}; do
   $srun python ncv_fold.py $path $i &
   # $srun python ncv_fold_with_correction.py $path $i &
done
wait
find $path -name preds.tsv | xargs tail -n +2 | grep -v "==>" | grep . | sort -k4 -n > $path/preds.tsv
$srun python collect_ncv_results.py $path

#sacct -j 15355615  --format=JobID,Start,End,Elapsed,REQCPUS,ALLOCTRES%30,Node

$srun python final_clock_cv.py $path