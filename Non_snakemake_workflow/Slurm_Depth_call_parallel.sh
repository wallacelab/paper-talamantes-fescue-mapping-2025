#!/bin/bash
#SBATCH -J SNP_Depth_P
#SBATCH -p batch
#SBATCH --ntasks=20
#SBATCH --mem 100gb
#SBATCH -t 160:00:00
#SBATCH --output=OutFiles/SNP_Depth_P.%j.out
#SBATCH -e OutFiles/SNP_Depth_P.%j.err
#SBATCH --mail-user=drt83172@uga.edu
#SBATCH --mail-type=FAIL,END

module load SAMtools/1.16.1-GCC-11.3.0
module load parallel/20221122-GCCcore-11.3.0
RefGenome="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Genome/tall_fescuev0.1.fa"
Bam_Files="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/bam_list.txt"
VCF_loc="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF"
bam_commands="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/bam_commands.txt"
# samtools depth --threads 4 -f $Bam_Files $VCF_loc/depths.bcf

parallel --jobs 3 < $bam_commands


