open(OUTA, ">$ARGV[0].slurm");

print OUTA "#!/bin/bash\n";
print OUTA "#SBATCH -J $ARGV[0]\n";
print OUTA "#SBATCH -p batch\n";
print OUTA "#SBATCH --ntasks=8\n";
print OUTA "#SBATCH --mem 40gb\n";
print OUTA "#SBATCH -t 75:00:00\n";
print OUTA "#SBATCH --output=OutFiles/$ARGV[0].%j.out\n";
print OUTA "#SBATCH -e OutFiles/$ARGV[0].%j.err\n";
print OUTA "#SBATCH --mail-type=NONE\n";
print OUTA "#SBATCH --mail-user=drt83172\@uga.edu\n";
print OUTA "#SBATCH --mail-type=FAIL\n";          

print OUTA "\n\n";

print OUTA "\npicard/2.27.5-Java-15\n";

print OUTA "java -jar \$EBROOTPICARD/picard.jar  MarkDuplicates \\
\tI=/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Parent_Bams/$ARGV[0]sorted.bam \\
\tO=/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/Parent_Bams/$ARGV[0]dupped.bam \\
\tM=/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Parent_Scripts/BWA_Scripts/Metrics/$ARGV[0]metrics.txt \\
\tREMOVE_DUPLICATES=true\n";




