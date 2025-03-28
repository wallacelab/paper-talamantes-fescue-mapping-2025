#!/bin/bash
#SBATCH -J Variant_Calling
#SBATCH -p batch
#SBATCH --ntasks=10
#SBATCH --mem 64gb
#SBATCH -t 48:00:00
#SBATCH --output=OutFiles/Variant_Calling.%j.out
#SBATCH -e OutFiles/Variant_Calling.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load Anaconda3/2022.10
source activate snakemake

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 32 -s snakefile

