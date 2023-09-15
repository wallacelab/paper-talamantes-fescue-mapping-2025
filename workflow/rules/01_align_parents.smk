rule bwa_map_parent:
    input:
        parent1=parentpath + "{parent_sample}R1_.fq.gz",
        parent2=parentpath + "{parent_sample}R2_.fq.gz",
        genome=config["genome"]
    output:
        bams2 = temp(mapped_reads_parents + "{parent_sample}.bam")
    threads:
        2
    conda:
        "../Conda_Envs/BWA.yaml"
    shell:
       "bwa mem {input.genome} {input.parent1} {input.parent2} -t {threads} > {output}"

rule samtools_sort_parent:
    input:
        bams2 = mapped_reads_parents + "{parent_sample}.bam"
    output:
        sorted2 = temp(mapped_reads_parents + "{parent_sample}sorted.bam")
    log:
        "logs/samtools/{parent_sample}_sorted.log"
    params:
        extra= "-m 16G",
    threads: 2
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule markduplicates_parents:
    input:
        bams = mapped_reads_parents + "{parent_sample}sorted.bam"
    output:
        bam = mapped_reads_parents + "{parent_sample}dupped.bam",
        metrics = picard_metrics + "{parent_sample}_metrics.txt"
    log:
        "logs/samtools/{parent_sample}_dupped.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"

# rule bcftools_mpileup_parent:
#     input:
#         alignments=mapped_reads_parents + "{parent_sample}dupped.bam",
#         ref=genome,
#         index=genome + ".fai",
#     output:
#         pileup= vcfs + "{parent_sample}.bcf",
#     params:
#         uncompressed_bcf=False,
#         extra="--max-depth 100 --min-BQ 10",
#     log:
#         "logs/bcftools_mpileup/{parent_sample}.log",
#     wrapper:
#         "v2.6.0/bio/bcftools/mpileup"
