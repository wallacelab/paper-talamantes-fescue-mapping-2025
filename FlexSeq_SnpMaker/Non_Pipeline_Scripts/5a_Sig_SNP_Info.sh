# This script is for after I identify the significant SNPs
# I found one significant SNP.
# 1. Using the SNP coordinates I removed it and 25,000 up and down stream nucleatides from the tall fescue genome
# 2. I make a blast data base and find genes in the extracted sequence


# This is the command used to extract the seqeunce
samtools faidx tall_fescue_pv1.1.fasta FaChr7G1:214286125-214336125 > SNP_FACHR7G1_214311125_sequence.fa
samtools faidx tall_fescue_pv1.1.fasta FaChr7G1:214301125-214321125 > SNP_FACHR7G1_214311125_sequence10k.fa


# This section will be all the blast commands I used 


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.2_Kyuss_2.0/GCF_019359855.2_Kyuss_2.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.2_Kyuss_2.0/GCF_019359855.2_Kyuss_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.2_Kyuss_2.0/GCF_019359855.2_Kyuss_2.0_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_protein.faa.gz 

# after files a gunzip you run this
makeblastdb -in GCF_019359855.2_Kyuss_2.0_genomic.fna -dbtype nucl -out lolium_perenne_Kyuss
blastn -query SNP_FACHR7G1_214311125_sequence.fa -db lolium_perenne_Kyuss -out results.out -outfmt 6


# Protein Blast
makeblastdb -in GCF_019359855.2_Kyuss_2.0_protein.faa -dbtype prot -out lp_Kyuss_protein_db
blastx -query SNP_FACHR7G1_214311125_sequence.fa -db lp_Kyuss_protein_db -outfmt "6 qseqid sseqid stitle pident length evalue bitscore" -evalue 1e-5 -out protein_hits.tsv
blastx -query SNP_FACHR7G1_214311125_sequence10k.fa -db lp_Kyuss_protein_db -outfmt "6 qseqid sseqid stitle pident length evalue bitscore" -evalue 1e-5 -out protein_hits_10k.tsv


# Rice Protein Blast
makeblastdb -in GCF_034140825.1_ASM3414082v1_protein.faa  -dbtype prot -out rice_db
blastx -query SNP_FACHR7G1_214311125_sequence.fa -db rice_db -outfmt "6 qseqid sseqid stitle pident length evalue bitscore" -evalue 1e-5 -out protein_hits_rice.tsv
blastx -query SNP_FACHR7G1_214311125_sequence10k.fa -db rice_db -outfmt "6 qseqid sseqid stitle pident length evalue bitscore" -evalue 1e-5 -out protein_hits_rice_10k.tsv



