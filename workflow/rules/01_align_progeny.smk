# aligns the files to the genome
rule bwa_map_prog:
    input:
        prog1=progpath + "{sample}R1.fastq.gz",
        prog2=progpath + "{sample}R2.fastq.gz",
        genome=config["genome"] 
    output:
        bams=temp(mapped_reads + "{sample}.bam")
    threads:
        2
    conda:
        "../Conda_Envs/BWA.yaml"
    shell:
       "bwa mem {input.genome} {input.prog1} {input.prog2} -t {threads} > {output}"



# Sorts the alignment file 
rule samtools_sort_prog:
    input:
        bams = mapped_reads + "{sample}.bam"   
    output:
        sorted = mapped_reads + "{sample}sorted.bam"
    log:
        "logs/samtools/{sample}_sorted.log",
    params:
        extra= "-m 2G",
    threads: 2
    wrapper:
        "v2.6.0/bio/samtools/sort"



# marks and removes the duplicates with picard
# rule markduplicates_prog:
#     input:
#         bams = mapped_reads + "{sample}sorted.bam"
#     output:
#         bam = mapped_reads + "{sample}dupped.bam",
#         metrics= picard_metrics + "{sample}_metrics.txt"
#     log:
#         "logs/samtools/{sample}_dupped.log",
#     params:
#         extra="--REMOVE_DUPLICATES true",
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v2.6.0/bio/picard/markduplicates"

