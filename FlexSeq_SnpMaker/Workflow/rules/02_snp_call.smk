# This rule should put all files into one huge bcf file
rule bcftools_mpileup_all:
    input:
        alignments = mapped_reads + "{sample}.sorted.bam",
        ref=genome,
        index=genome + ".fai",
    output:
        pileup= mpileup + "{sample}_pileup.bcf",
    params:
        uncompressed_bcf=False,
        extra="--max-depth 100 --min-BQ 10",
    log:
        "logs{sample}_mpileup_call.log",
    wrapper:
        "v2.6.0/bio/bcftools/mpileup"


rule bcftools_call:
    input:
        pileup = mpileup + "{sample}_pileup.bcf"
    output:
        vcf = vcfs + "{sample}.vcf"
    log:
        "logs/bcftools_call{sample}.log"
    conda:
        "../Conda_Envs/bcftools.yaml"
    shell:
        """
        bcftools call -mv {input.pileup} -o {output}
        """








