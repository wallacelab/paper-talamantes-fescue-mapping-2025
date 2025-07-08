# Tall fescue mapping & QTL analysis
*Repo of scripts for QTL mapping in Tall Fescue, by Darrian Talamantes*

This will filter reads, align my files to the tall fescue genome, find SNPs, create a genetic map 

Programs and what they do:

These SLURM_GENERATOR scripts are made to work with parent files right now.
- SLURM_GENERATOR_BWA.pl: runs BWA on all files
- SLURM_GENERATOR_Sorting_BAMS.pl: Sorts the bam files
- SLURM_GENERATOR_samtools_idx: indexes bam files
- SLURM_GENERATOR_Remove_Dupes.pl: removes duplicates from bam files
- SLURM_GENERATOR_samtools_idx.pl: indexes the files after removing duplicates
- SLURM_GENERATOR_parent_bam_filter.pl: Removes all regions not found in progeny files


The workflow folder contains a snakemake pipline used to process all the progeny files locally.


R files to know:
- phenotype_analysis.R: This file takes the CT data and processes it. It also process alkaloid data. It will then create risidual data. Filtered CT refers to CT values after accounting for efficiency of each run.
- Adds_paternal_data_to_files.R: This file adds the metadata to the residual phenotype file.
- Allel_frequincy_evaluation.R: Evaluates how normal the allel frequncies are.
- VCF_Analysis_Manhattan_Plots.R: Analyses depth and takes GWAS results to make manhattan plots.


Important Data Files:
- Variants_all_ez_names.bcf: Variants from all files in one place with names that dont include directories.
- residual_data_father_included_2.txt: residual data with all the metadata
- Alkaloid_Data_Combined.csv: This has the alkaloid data with its metadata
- Phenotype_Data_Delta_CT.txt: This has the biomass data delta CT only. It is created by phenotype analysis if we want the ratio of fescue/epi in the future.


