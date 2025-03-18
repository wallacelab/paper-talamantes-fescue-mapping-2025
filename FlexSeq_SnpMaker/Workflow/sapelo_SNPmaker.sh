#!/bin/bash
#SBATCH -J FlesSeqSNPMaker
#SBATCH -p batch
#SBATCH --ntasks=32
#SBATCH --mem 120gb
#SBATCH -t 160:00:00
#SBATCH --output=OutFiles/FlesSeqSNPMaker.%j.out
#SBATCH -e OutFiles/FlesSeqSNPMaker.%j.err
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

