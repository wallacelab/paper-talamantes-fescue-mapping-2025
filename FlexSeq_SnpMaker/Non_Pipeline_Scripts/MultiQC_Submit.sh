#!/bin/bash
#SBATCH -J MultiQC
#SBATCH -p highmem_p
#SBATCH --ntasks=32
#SBATCH --mem 500gb
#SBATCH -t 140:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Outfiles/MultiQC.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Outfiles/MultiQC.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load MultiQC/1.14-foss-2022a

multiqc .