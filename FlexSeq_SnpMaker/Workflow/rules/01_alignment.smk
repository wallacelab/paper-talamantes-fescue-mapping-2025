# aligns the files to the genome
rule bwa_map:
    input:
        sample=samples_path + "{sample}fasta",
        genome=config["genome"],
        idx = multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa") 
    output:
        bams=temp(mapped_reads + "{sample}.bam")
    threads:
        2
    conda:
        "../Conda_Envs/BWA.yaml"
    shell:
       "bwa mem {input.genome} {input.sample} -t {threads} > {output}"



# Sorts the alignment file 
rule samtools_sort:
    input:
        bams = mapped_reads + "{sample}.bam"   
    output:
        sorted = mapped_reads + "{sample}sorted.bam"
    log:
        "logs/samtools{sample}_sorted.log",
    params:
        extra= "-m 2G",
    threads: 2
    wrapper:
        "v2.6.0/bio/samtools/sort"

