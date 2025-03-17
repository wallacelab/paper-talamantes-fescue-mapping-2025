rule align_reads:
    input:
        ref = config["genome"],  # Reference genome
        refidx1 = multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        refidx2 = genome + ".fai",
        fq1 = samples_path + "/{sample}_R1.fastq.gz",
        fq2 = samples_path + "/{sample}_R2.fastq.gz"
    output:
        bam = temp(mapped_reads + "/{sample}.bam")
    log:
        "logs/{sample}_align.log"
    threads: 16
    conda:
        "../Conda_Envs/BWA.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fq1} {input.fq2} -o {output.bam} 2>> {log}
        """

# Sorts the alignment file 
rule samtools_sort:
    input:
        bams = mapped_reads + "/{sample}.bam"   
    output:
        sorted = mapped_reads + "/{sample}_sorted.bam"
    log:
        "logs/samtools{sample}_sorted.log",
    params:
        extra= "-m 2G",
    threads: 16
    wrapper:
        "v2.6.0/bio/samtools/sort"

