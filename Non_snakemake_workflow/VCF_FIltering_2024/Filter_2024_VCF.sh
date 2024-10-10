#Objective: This will filter the snps from Freedom markers 2024 data.
list_314x312=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/314x312_2024.txt
list_314x310=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/314x310_2024.txt

vcf_file=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.all.snps.vcf.gz
vcf_mcr50=/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.MCR50.snps.vcf.gz

Major=0.85
missing=43
depth=696
# Filters I am using
echo "filtering MAF at $MAF"
echo "filtering missing at $missing"
echo "filtering depth at $depth" 
bcftools view -s `cat $list_314x312 | tr '\n' ',' | sed 's/,$//'` $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf
bcftools view -s `cat $list_314x312 | tr '\n' ',' | sed 's/,$//'` $vcf_mcr50 --force-samples -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x312.vcf


echo "filtering MAF at $MAF"
echo "filtering missing at $missing"
echo "filtering depth at $depth" 
bcftools view -s `cat $list_314x310 | tr '\n' ',' | sed 's/,$//'` $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf
bcftools view -s `cat $list_314x310 | tr '\n' ',' | sed 's/,$//'` $vcf_mcr50 --force-samples -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x310.vcf

echo "starting $list_314x312"
# bcftools view -S $list_314x312 $vcf_file -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf --maf 0.15 --max-missing 0.5 --minDP 8 --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered

echo "starting $list_314x310"
vcftools --vcf /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf --maf 0.15 --max-missing 0.5 --minDP 8 --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered

echo "starting $vcf_mcr50"
vcftools --vcf $vcf_mcr50 --maf 0.15 --max-missing 0.5 --minDP 8 --recode --out /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_Extra_filtered

rm /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310.vcf
rm /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf









################ Below is kinda useless tbh ########################################


# echo "bcftools view -i INFO/AF[0] > $Major && INFO/DP > $depth && F_MISSING && INFO/NS >= $missing /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered.vcf"
# bcftools view -i "INFO/AF[0] <= $Major && INFO/DP >= $depth && F_MISSING && INFO/NS >= $missing" /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312.vcf -o /home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered.vcf


########## Examples of filtering with bcftools
# To filter by depth we have 87 samples in cross 312x314. For avg depth of 10 we would use 870
# -i 'INFO/DP > 870'

# Major allel frequency filter 
# 'INFO/AF[0] > 0.85'

# Get rid of genotypes missing more than 50% 43/87
# -i 'INFO/NS >= 43' 

# Keep only bi allilic sites
# -m2 -M2

# bcftools view -s $list_314x310



