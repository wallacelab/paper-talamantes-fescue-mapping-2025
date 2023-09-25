open(OUTA, ">$ARGV[0].slurm");

print OUTA "#!/bin/bash\n";
print OUTA "#SBATCH -J $ARGV[0]\n";
print OUTA "#SBATCH -p batch\n";
print OUTA "#SBATCH --ntasks=8\n";
print OUTA "#SBATCH --mem 40gb\n";
print OUTA "#SBATCH -t 75:00:00\n";
print OUTA "#SBATCH --output=$ARGV[0].%j.out\n";
print OUTA "#SBATCH -e $ARGV[0].%j.err\n";
print OUTA "#SBATCH --mail-type=NONE\n";
print OUTA '#SBATCH --mail-user drt83172@uga.edu';
print OUTA "#SBATCH --mail-type=FAIL";          
print OUTA "\n\n";

print OUTA "\nmodule load BWA/0.7.17-GCCcore-11.2.0\n";
print OUTA "bwa mem /scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Genome/tall_fescuev0.1.fa /scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Parent_Reads/$ARGV[0]_R1_.fq.gz /scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Parent_Reads/$ARGV[0]_R2_.fq.gz -t 8 > /scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Parent_Bams/$ARGV[0].bam\n";
