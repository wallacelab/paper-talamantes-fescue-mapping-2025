#Objective: This will filter the snps from Freedom markers 2024 data.
list_314x312=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/314x312_2024.txt
list_314x310=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/314x310_2024.txt

vcf_file=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.all.snps.vcf.gz
vcf_mcr50=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.MCR50.snps.vcf.gz

MAF=0.05
missing=.5
depth=8
maxdepth=30

#Output files
$all_output=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered.recode.vcf

###################################
# This extracts populations from larger files
###################################

echo
echo "creating subseted data for 314x312"
echo 
bcftools view -s `cat $list_314x312 | tr '\n' ',' | sed 's/,$//'` $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf
bcftools view -s `cat $list_314x312 | tr '\n' ',' | sed 's/,$//'` $vcf_mcr50 --force-samples -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x312.vcf


echo
echo "creating subseted data for 314x310"
echo 
bcftools view -s `cat $list_314x310 | tr '\n' ',' | sed 's/,$//'` $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf
bcftools view -s `cat $list_314x310 | tr '\n' ',' | sed 's/,$//'` $vcf_mcr50 --force-samples -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x310.vcf

##################################
# Filters extracted data and the all snps file with my filters
##################################
echo 
echo "############################################################################"
echo "filtering MAF at $MAF"
echo "filtering missing at $missing"
echo "filtering depth at $depth and max of $maxdepth"
echo "############################################################################"
echo 

echo "Filtering $list_314x312"
# bcftools view -S $list_314x312 $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf --maf $MAF --max-missing $missing --minDP $depth --maxDP $maxdepth --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered

echo "Filtering $list_314x310"
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf --maf $MAF --max-missing $missing --minDP $depth --maxDP $maxdepth --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered

echo "Filtering $vcf_file"
vcftools --gzvcf $vcf_file --maf $MAF --max-missing $missing --minDP $depth --maxDP $maxdepth --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered

rm /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf
rm /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf

##################################
# Getting rid of taxa with too much missing
##################################
echo 
echo #################################################
echo "Getting rid of taxa with too much missing"
echo #################################################
echo

# all_snps_filtered
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered.recode.vcf --missing-indv --out missing_data_report
awk '$5 > 0.9 {print $1}' missing_data_report.imiss > taxa_to_remove.txt
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered.recode.vcf --remove taxa_to_remove.txt --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2

# all_snps_filtered 314x310
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered.recode.vcf --missing-indv --out missing_data_report
awk '$5 > 0.9 {print $1}' missing_data_report.imiss > taxa_to_remove.txt
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered.recode.vcf --remove taxa_to_remove.txt --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered_2

# all_snps_filtered 314x310
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered.recode.vcf --missing-indv --out missing_data_report
awk '$5 > 0.9 {print $1}' missing_data_report.imiss > taxa_to_remove.txt
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered.recode.vcf --remove taxa_to_remove.txt --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered_2


#######
# Getting all the files needed for vcf file stat analysis
#######

# Depth files
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf --site-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_site_depth

vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered_2.recode.vcf --site-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_site_depth

vcftools --gzvcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.MCR50.snps.vcf.gz --site-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_site_depth

# Genotype Depth file
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf --geno-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_depth

vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered_2.recode.vcf --geno-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_depth

vcftools --gzvcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.MCR50.snps.vcf.gz --geno-depth --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_depth


#########
# List of final output vcf files for my filters
#########
# /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf
# /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered_2.recode.vcf
# /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered_2.recode.vcf



######
# Below is me running the R program that makes summary stats for my VCF files
######

#./1c_PlotGenoSummary.r -d /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_site_depth.ldepth -g /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_depth.gdepth -t /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_TaxaSummary.txt -s /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_SiteSummary.txt -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_filtered_v2.png

# ./1c_PlotGenoSummary.r -d /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_site_depth.ldepth -g /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_depth.gdepth -t /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_TaxaSummary.txt -s /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_SiteSummary.txt -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/MCR50_stats.png


#./1c_PlotGenoSummary.r -d /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_site_depth.ldepth -g /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_depth.gdepth -t /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_TaxaSummary.txt -s /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered_SiteSummary.txt -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/files_for_site_analysis/all_snps_314x312_filtered.png




