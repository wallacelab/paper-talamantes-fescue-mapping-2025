echo 
2 genome
3 threads



for $file in $(ls $1)
do 
    arrIN=(${file//./ }) # makes the variable into an array that I sepetate by "."
    prog=${arrIN[0]}
    bwa mem $2 $file -t $3 > mapped_reads/$prog.bam
done 




# sync_FLEX-SEQ_RAW_UGA_149001_P001_WA01_RAPiD-Genomics_F165_UGA_149001_P001_WA01_i5-54_i7-59_S1_L002_R1_001.fastq.gz
