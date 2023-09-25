#!/bin/bash
#SBATCH -J SNP_CALL_315x320
#SBATCH -p batch
#SBATCH --ntasks=32
#SBATCH --mem 120gb
#SBATCH -t 80:00:00
#SBATCH --output=OutFiles/SNP_CALL.%j.out
#SBATCH -e OutFiles/SNP_CALL.%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=drt83172@uga.edu
#SBATCH --mail-type=FAIL

module load BCFtools/1.15.1-GCC-11.3.0
RefGenome="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Genome/tall_fescuev0.1.fa"
Bam_Files="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/Full_sibs/315x320.txt"
VCF_loc="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF"

bcftools mpileup -f $RefGenome -b $Bam_Files | bcftools call -mv -Ob -o $VCF_loc/315x320_Variants.bcf
