#!/bin/bash
#SBATCH -J filter_312x314
#SBATCH -p batch
#SBATCH --ntasks=8
#SBATCH --mem 20gb
#SBATCH -t 24:00:00
#SBATCH --output=OutFiles/SNP_CALL.%j.out
#SBATCH -e OutFiles/SNP_CALL.%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=drt83172@uga.edu
#SBATCH --mail-type=FAIL

#Module loading
module load VCFtools/0.1.16-GCC-11.2.0
module load BCFtools/1.15.1-GCC-11.3.0

vcf_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF"

vcftools --bcf $vcf_loc/312x314_Variants.bcf --maf .05 --max-maf .95 --max-alleles 2 --exclude-uncalled  --out $vcf_loc/312x314_filtered

# --max-af 95
# --min-af .05 --max-alleles 2
# --include 
