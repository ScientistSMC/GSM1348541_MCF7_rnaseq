#!/bin/bash
#SBATCH --job-name read_counts ##name that will show up in the queue
#SBATCH --output read_counts.o%j ##filename of the output; the "%j" will append the jobID to the end of the name making the output files unique despite the sane job name; default is slurm-[jobID].out
#SBATCH --partition normal ##the partition to run in [options: normal, gpu, debug]; default = normal
#SBATCH	--nodes 1 ##number of nodes to use; default = 1
#SBATCH --ntasks 4 ##number of tasks (analyses) to run; default = 1
#SBATCH --cpus-per-task 8 ##the number of threads the code will use; default = 1
#SBATCH --mem 100000
#SBATCH	--time 0-12:00:00 ##time for analysis (day-hour:min:sec) -- Max walltime will vary by partition; time formats are: "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#SBATCH --mail-user samcarli@nmsu.edu   ##your email address
#SBATCH --mail-type BEGIN ##slurm will email you when your job starts
#SBATCH --mail-type END ##slurm will email you when your job ends
#SBATCH --mail-type FAIL ##slurm will email you when your job fails
#SBATCH --get-user-env ##passes along environmental settings 

module load anaconda

conda activate R

Rscript 2a.compute_summarizeOverlaps.ebg.strict.R

