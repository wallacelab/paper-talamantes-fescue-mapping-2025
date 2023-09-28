
vcf_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF"

bcftools view --min-af .05 --min-alleles 2 --exclude-uncalled $vcf_loc/Variants_progeny.bcf > $vcf_loc/Variants_progeny.vcf

# --max-af 95
# --min-af .05 --min-alleles 2
# --include 
