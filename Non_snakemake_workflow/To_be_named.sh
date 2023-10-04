#!/bin/bash

# Author: Darrian Talamantes
# Purpose: This will filter vcf files in ways where I can then import the data into R for analysis and graphs
# The R script that goes with this is named: VCF_Analysis.R

VCF_file=/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/Variants_all.vcf


# To pull depth from VCF.
cut -f 1,2,8 $VCF_file > $VCF_file.extract
sed 's/\;.*//' $VCF_file.extract > $VCF_file.extract.1
sed 's/\DP=//' $VCF_file.extract.1 > $VCF_file.extract.counts
