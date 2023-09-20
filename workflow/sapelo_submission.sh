#!/bin/bash
#SBATCH --job-name=Variant_Calling
#SBATCH --partition=batch 
#SBATCH --nodes=1 
#SBATCH --ntasks=32
#SBATCH --time=120:00:00
#SBATCH --mem=110gb
#SBATCH --mail-user=drt83172@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=OutFiles/Variant_Calling.%j.out
#SBATCH -e OutFiles/Variant_Calling.%j.err

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 20
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load Anaconda3/2022.10
source activate snakemake

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 32 -s snakefile

