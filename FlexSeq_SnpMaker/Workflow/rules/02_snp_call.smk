# This rule should put all files into one huge bcf file
rule call_snps:
    input:
        sorted = expand(mapped_reads + "/{sample}_sorted.bam", sample=SAMPLES),
        ref=genome,
        index=genome + ".fai"
    output:
        vcf= vcfs + "/all_Flexseq.vcf.gz"
    log:
        "logs/call_snps.log"
    threads: 16
    conda:
        "../Conda_Envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} {input.sorted} --threads {threads} | bcftools call -mv -Oz -o {output.vcf} --threads {threads} 2>> {log}
        bcftools index {output.vcf}
        """






