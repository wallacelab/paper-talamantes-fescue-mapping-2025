#!/bin/bash
#SBATCH -J FlexSeqSNPMaker
#SBATCH -p highmem_p
#SBATCH --ntasks=16
#SBATCH --mem 300gb
#SBATCH -t 160:00:00
#SBATCH --output=OutFiles/FlexSeqSNPMaker.%j.out
#SBATCH -e OutFiles/FlexSeqSNPMaker.%j.err
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

snakemake --use-conda --cores 16 -s snakefile

