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
    threads: 8
    conda:
        "../Conda_Envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} {input.sorted} | bcftools call -mv -Oz -o {output.vcf}
        if [[ $? -ne 0 ]]; then echo "bcftools call failed" >&2; exit 1; fi
        bcftools index {output.vcf}

        """






